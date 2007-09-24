#ifndef CRANKNICHOLSON_H
#define CRANKNICHOLSON_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/cartesianrepresentation.h"
#include "../potential/examplepotentials.h"
#include "../utility/fortran.h"

extern "C"
{
	void FORTRAN_NAME(zgtsv)( int *N, int *nhrs, cplx *diagLower, cplx* diag, cplx* diagUpper, cplx* B, int* LDB, int *info );
};

/*
 * Propagates the potential and kinetic energy using the crank nicholson
 * scheme
 *
 * It is implemented as a potential evaluator. That way we can specify a potential 
 * which it will propagate without the need of a split step propagation."
 */


/*
 * DynamicPotentialClass provides the method GetPotential(pos,t, dt)
 */


template<class DynamicPotentialClass, int Rank>
class CrankNicholsonEvaluator : public DynamicPotentialClass
{

public:
	typedef blitz::Array<cplx, 1> DataArray;

	int RadialRank;
	double Mass;

	CrankNicholsonEvaluator() : RadialRank(0), Mass(1.0) {}

	/* 
	 * Reads radial rank (the rank which this CN-evaluator will propagate)
	 * and mass from the configuration section of the potential
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		DynamicPotentialClass::ApplyConfigSection(config);
		if (config.HasValue("radial_rank"))
		{
			config.Get("radial_rank", RadialRank);
		}

		if (config.HasValue("mass"))
		{
			config.Get("mass", Mass);
		}
	}


	/** Propagates the wavefunction a part of one timestep this function should be called
	 * with parity = 0 and parity = 1 in order to complete a timestep. Ghost points must be updated
	 * between calls with different parity
	 **/
	void UpdateWavefunction(Wavefunction<Rank> &psi, double t, cplx dt)
	{
		//Set up t and dt for the potential
		this->TimeStep = dt;
		this->CurTime = t;
		
		//Get dx
		CartesianRepresentation<1>* repr = this->GetRepresentation(psi);
		double dx = repr->GetRange(0).Dx;
		double startx = repr->GetRange(0).Min;
		int gridSize = repr->GetRange(0).Count;

		//C-N scheme for schrod: (i - dt/2 A)U^(n+1) = (i + dt/2 A)U^(n+1)
		//A = - 1 / (2 m) D^2 + V
		//i is the imaginary unit
		//D^2 is the 2. order central finite difference approx to the 2. derivative
		//V is the potential
	
		//the fwd matrix is (i + dt/2 A) = (i - dt / (4 m) D^2) 
		//the bwd matrix is (i - dt/2 A) = (i + dt / (4 m) D^2)
		
		//Forward subdiagonal
		DataArray fwdDiagonal(gridSize);
		DataArray fwdUpperDiagonal(gridSize-1); 
		DataArray fwdLowerDiagonal(gridSize-1); 
		fwdUpperDiagonal = - dt / (4.0 * Mass * dx * dx);
		fwdLowerDiagonal = - dt / (4.0 * Mass * dx * dx);

		//Backward subdiagonal
		DataArray bwdDiagonal(gridSize);
		DataArray bwdUpperDiagonal(gridSize-1);
		DataArray bwdLowerDiagonal(gridSize-1);
		bwdUpperDiagonal = dt / (4.0 * Mass * dx * dx);
		bwdLowerDiagonal = dt / (4.0 * Mass * dx * dx);

		//Setup diagonals
		blitz::TinyVector<double, 1> pos(startx);
		const cplx imaginaryUnit(0.0, 1.0);
		for (int i=0; i<gridSize; i++)
		{
			cplx pot = this->GetPotentialValue(pos);
			pos(0) += dx;
			cplx H = (2.0 / ( 2.0 * Mass * dx * dx)) + pot;
			fwdDiagonal(i) = imaginaryUnit + (dt / 2.0) * H;
			bwdDiagonal(i) = imaginaryUnit - (dt / 2.0) * H;
		}

		//propagate the wavefunction
		DataArray psiData(psi.GetData());
		blitz::Array<cplx, 1> tempData(gridSize);

		tempData = psiData;
		SolveForward(fwdDiagonal, fwdUpperDiagonal, fwdLowerDiagonal, tempData, psiData);
		SolveBackward(bwdDiagonal, bwdUpperDiagonal, bwdLowerDiagonal, psiData);
	}


	/*
	 * Applies the differentiation operator + potential operator to the wavefunction
	 * useful for calculating expecation values and possibly finding eigenstates
	 */
	void MultiplyOperator(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi, double t, const cplx &dt)
	{
		//Set up t and dt for the potential
		this->TimeStep = dt;
		this->CurTime = t;

		//Get dx
		CartesianRepresentation<1>* repr = this->GetRepresentation(srcPsi);
		double dx = repr->GetRange(0).Dx;
		double startx = repr->GetRange(0).Min;
		int gridSize = repr->GetRange(0).Count;

		blitz::Array<cplx, 1> srcData = srcPsi.GetData();
		blitz::Array<cplx, 1> dstData = dstPsi.GetData();

		blitz::TinyVector<double, 1> pos(startx);
		cplx prevData = 0;
		//By taking curData=0, and nextData=data(0) we eliminate the first boundary condition
		cplx curData = 0; 
		cplx nextData = srcData(0);
		double scaling = - 1.0 / (2.0 * Mass * dx * dx);
		for (int i=0; i<gridSize-1; i++)
		{
			//1D 2. order central difference
			prevData = curData;
			curData = nextData;
			nextData = srcData(i+1);
			cplx finiteDiff = (prevData - 2.0 * curData + nextData) * scaling;

			//Potential
			cplx pot = this->GetPotentialValue(pos) * curData; 
			pos(0) += dx;

			dstData(i) += pot + finiteDiff;
		}
		//The last grid point
		dstData(gridSize-1) += (curData - 2.0 * nextData) * scaling + this->GetPotentialValue(pos) * curData;
	}


private:
	CartesianRepresentation<1>* GetRepresentation(Wavefunction<Rank> &psi)
	{
		return (CartesianRepresentation<1>*)&psi.GetRepresentation();
	}

	void SolveForward(DataArray &diag, DataArray &upperDiag, DataArray &lowerDiag, DataArray &in, DataArray &out)
	{
		int N = diag.extent(0);
		for (int i=1; i<N-1; i++)
		{
			out(i) = lowerDiag(i-1) * in(i-1) 
                      + diag(i)      * in(i) 
                      + upperDiag(i) * in(i+1);
		}

		out(0) = diag(0) * in(0) + upperDiag(0) * in(1);
		out(N-1) = lowerDiag(N-2) * in(N-2) + diag(N-1) * in(N-1);
	}

	void SolveBackward(DataArray diag, DataArray upperDiag, DataArray lowerDiag, DataArray data)
	{
		int N = diag.extent(0);
		int rhsCount = 1;
		int leadingB = N;
		int info = 0;
		FORTRAN_NAME(zgtsv)( &N, &rhsCount, lowerDiag.data(), diag.data(), upperDiag.data(), data.data(), &leadingB, &info);

		if (info != 0)
		{
			cout << "Error solving backward triangular system " << info << endl;
		}
	}
};


#endif


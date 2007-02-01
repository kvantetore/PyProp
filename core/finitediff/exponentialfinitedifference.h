#ifndef EXPONENTIALFINITEDIFFERENCE_H
#define EXPONENTIALFINITEDIFFERENCE_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/cartesianrepresentation.h"
#include "../potential/examplepotentials.h"

/*
 * - Propagates one timestep of a potential and kinetic energy.
 * - Uses finite difference to approximate potential energy operator
 *   this gives a tri diagonal matrix.
 * - Uses split step, to split the tri diagonal matrix into two block 
 *   diagonal matrices (of parity 0 and 1)  
 * - The exponential of a tri diagonal matrix can easily be found in
 *   O(n) operations, by diagonalizing each 2x2 block seperatly
 */


/*
 * PotentialEvaluator provides the method GetPotential(pos,t, dt)
 */


template<class DynamicPotentialClass, int Rank>
class ExponentialFiniteDifferenceEvaluator : public DynamicPotentialClass
{
public:

	/** Propagates the wavefunction a part of one timestep this function should be called
	 * with parity = 0 and parity = 1 in order to complete a timestep. Ghost points must be updated
	 * between calls with different parity
	 **/
	void UpdateWavefunction(Wavefunction<Rank> &psi, double t, cplx dt, int parity)
	{
		//Get dx
		CartesianRepresentation<1>* repr = (CartesianRepresentation<1>*)&psi.GetRepresentation();
		double dx = repr->GetRange(0).Dx;
		double startx = repr->GetRange(0).Min;

		//Set up t and dt for the potential
		this->TimeStep = dt;
		this->CurTime = t;
	
		//Get this procs position in the global scheme of things...
		DistributedModel<1>* distr = & psi.GetRepresentation().GetDistributedModel();
		int globalStartIndex = distr->GetGlobalStartIndex(psi,0);
	
		//If we do not start on the beginning of a block
		int startIndex = 0;
		if (globalStartIndex % 2 != parity)
		{
			blitz::TinyVector<double, Rank> pos;
			pos(0) = startx;
			double pot2 = GetPotentialValue(pos);
			
			if (distr->IsFirstProc())
			{
				//Apply diagonal
				EXPFD_UpdateDiagonal(pot2, psi.Data(0), dx, dt);
			}
			else
			{
				//Apply left ghost block
				pos(0) -= dx;
				//double pot2 = GetPotentialValue(pos);
			
				throw std::runtime_error("Not implemented ghost points yet");
			}
			
			startIndex++;
			startx += dx;
		}
	
		//If we do not end on the end of a block
		int lastIndex = psi.Data.extent(0) - 1;
		if ((globalStartIndex + lastIndex) % 2 == parity)
		{
			blitz::TinyVector<double, Rank> pos;
			pos(0) = startx + dx * lastIndex;
			double pot1 = GetPotentialValue(pos);
			
			if (distr->IsLastProc())
			{
				//Apply diagonal
				EXPFD_UpdateDiagonal(pot1, psi.Data(lastIndex), dx, dt);
			}
			else
			{
				pos(0) += dx;
				//double pot2 = GetPotentialValue(pos);
			
				//Apply right ghost block
				throw std::runtime_error("Ghost points not implemented yet");
			}
			
			lastIndex--;
		}
	

		
		//All inner points
		blitz::TinyVector<double, Rank> pos;
		pos(0) = startx;
		for (int i=startIndex; i<=lastIndex; i+=2)
		{
			//Apply block
			double pot1 = GetPotentialValue(pos);
			pos(0) += dx;
			double pot2 = GetPotentialValue(pos);
			pos(0) += dx;

			EXPFD_UpdateBlock(pot1, pot2, psi.Data(i), psi.Data(i+1), dx, dt);
		}
	}
};

/* Local function definitions. */
double inline EXPFD_GetDiagonalValue(const double &pot, const double &dx);
void inline EXPFD_UpdateEigenvalues(const double &a, const double &b, const double &c, double &l1, double &l2);
void inline EXPFD_UpdateEigenvectors(const double &a, const double &b, const double &c, double &l1, double &l2, blitz::TinyVector<double, 2> &v1, blitz::TinyVector<double, 2> &v2);
void inline EXPFD_UpdateBlock(const double &pot1, const double &pot2, cplx &psi1, cplx &psi2, const double &dx, const cplx &dt);
void inline EXPFD_UpdateDiagonal(const double &pot, cplx &psi, const double &dx, const cplx &dt);

/** Gets the eigenvalues for the matrix ([[a, b], [b, c]] */
void inline EXPFD_UpdateEigenvalues(const double &a, const double &b, const double &c, double &l1, double &l2)
{
	double subval = sqrt( sqr(a+c) - 4 * (a*c - sqr(b)));
	l1 = 0.5 * (a + c + subval);
	l2 = 0.5 * (a + c - subval);
}

/** Gets the eigenvectors and eigenvalues for the matrix [[a, b], [b, c]] **/
void inline EXPFD_UpdateEigenvectors(const double &a, const double &b, const double &c, double &l1, double &l2, blitz::TinyVector<double, 2> &v1, blitz::TinyVector<double, 2> &v2)
{
	EXPFD_UpdateEigenvalues(a,b,c,l1,l2);

	double direction = b / (l1 - a);
	double lengthDiv = sqrt(1.0 / (1 + sqr(direction)));
	v1(0) = direction * lengthDiv;
	v1(1) = lengthDiv;

	direction = b / (l2 - a);
	lengthDiv = sqrt(1.0 / (1 + sqr(direction)));
	v2(0) = direction * lengthDiv;
	v2(1) = lengthDiv;
}

/** Get the diagoal value of the matrix **/
double inline EXPFD_GetDiagonalValue(const double &pot, const double &dx)
{
	return - 0.5 * (pot * sqr(dx) + 1.0);
}

/*
 * Given the 2x2 matrix B = [[a,b],[b,c]]
 * Where a, b and c are given by the potential and dx,
 * We will compute psi = psi * exp(i * B * dt)
 */
void inline EXPFD_UpdateBlock(const double &pot1, const double &pot2, cplx &psi1, cplx &psi2, const double &dx, const cplx &dt)
{
	double a = EXPFD_GetDiagonalValue(pot1, dx);
	double b = 0.5;
	double c = EXPFD_GetDiagonalValue(pot2, dx);
	
	double l1, l2;
	blitz::TinyVector<double, 2> v1, v2;
	EXPFD_UpdateEigenvectors(a, b, c, l1, l2, v1, v2);

	// psi(t+dt) = exp(A) psi = U exp(D) U* psi 
	// U  = [[v1(0) , v2(0) ], [v1(1) , v2(1) ]]
	// U* = [[v1(0)*, v1(1)*], [v2(0)*, v2(1)*]]
	// D  = [[l1    , 0     ], [0     , l2    ]] * I * (dt/dx**2)
	
	//t = U* psi
	cplx t1 = v1(0) * psi1 + v1(1) * psi2;
	cplx t2 = v2(0) * psi1 + v2(1) * psi2;

	//t *= exp(D)
	t1 *= exp( I * (dt / sqr(dx)) * l1);
	t2 *= exp( I * (dt / sqr(dx)) * l2);

	//psi(t+dt) = U t
	psi1 = v1(0) * t1 + v2(0) * t2;
	psi2 = v1(1) * t1 + v2(1) * t2;
}

/* Updates the diagonal elements (first/last) */
void inline EXPFD_UpdateDiagonal(const double &pot, cplx &psi, const double &dx, const cplx &dt)
{
	psi = psi * exp(I * (dt / sqr(dx)) * EXPFD_GetDiagonalValue(pot, dx));
}


#endif


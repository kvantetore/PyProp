#ifndef PAMPWRAPPER_H
#define PAMPWRAPPER_H

#include <core/common.h>
#include <core/wavefunction.h>
#include <core/utility/boostpythonhack.h>

#include "pamp/pamp.h"
#include "../pypropfunctor.h"

namespace krylov
{

template<int Rank>
class PampWrapper : public PypropKrylovWrapper
{
public:
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;

private:
	pamp::pAMP<cplx> Propagator;

	//Temporary member variables
	double CurTime;
	cplx TimeStep;
	typename Wavefunction<Rank>::Ptr Psi;
	typename Wavefunction<Rank>::Ptr TempPsi;
	object Callback;

	//test
	//blitz::Array<cplx, 2> M;

public:
	void ApplyConfigSection(const ConfigSection &config);
	void Setup(const typename Wavefunction<Rank>::Ptr psi);
	void AdvanceStep(object callback, typename Wavefunction<Rank>::Ptr psi, typename Wavefunction<Rank>::Ptr tempPsi, cplx dt, double t, bool usePypropIntegration);

	void ApplyOperator(blitz::Array<cplx, 1> &input, blitz::Array<cplx, 1> &output);

	/*
	 * Returns the _converged_ eigenvectors, a N by M matrix, where
	 * N == GetEigenvalueCount(), and M is the size of the wavefunction
	 */
	double EstimateMemoryUsage(int matrixSize, int basisSize)
	{
		return Propagator.EstimateMemoryUsage(matrixSize, basisSize);
	}

	/* 
	 * Returns the nuber of matrix-vector operations performed
	 */
	int GetOperatorCount()
	{
		return Propagator.GetOperatorCount();
	}

	/*
	 * Returns the number of re-orthogonalization steps performed
	 */
	int GetOrthogonalizationCount()
	{
		return Propagator.GetOrthogonalizationCount();
	}

	void PrintStatistics()
	{
		Propagator.PrintStatistics();
	}

	void ResetStatistics()
	{
		Propagator.ResetStatistics();
	}

	void TestPade()
	{
		typedef typename blitz::linalg::LAPACK<cplx>::MatrixType MatrixType;
		typedef typename blitz::linalg::LAPACK<cplx>::VectorType VectorType;

		MatrixType matrix(2,2);
		MatrixType exp(matrix.shape());
		matrix(0,0) = -49;
		matrix(1,0) = 24;
		matrix(0,1) = 24;
		matrix(1,1) = 31;

		Propagator.PadeExponential(matrix, exp, 15);
		cout << "A = " << matrix << endl;
		cout << "exp(A) = " << exp << endl;
		
		matrix *= cplx(0,1);
		Propagator.PadeExponential(matrix, exp, 15);
		cout << "exp(iA) = " << exp << endl;

		Propagator.ScalingPadeExponential(matrix, exp, -1, -1);
		cout << "exp(iA) = " << exp << endl;
	}

	double GetResidualNorm()
	{
		return Propagator.GetResidualNorm();
	}

};

} // Namespace

#endif


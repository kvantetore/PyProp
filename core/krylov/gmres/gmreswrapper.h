#ifndef PAMPWRAPPER_H
#define PAMPWRAPPER_H

#include <core/common.h>
#include <core/wavefunction.h>
#include <core/utility/boostpythonhack.h>

#include "gmres.h"
#include "../pypropfunctor.h"

namespace krylov
{

template<int Rank>
class GmresWrapper : public PypropKrylovWrapper
{
public:
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;

private:
	gmres::GMRES<cplx> Solver;

	//Temporary member variables
	typename Wavefunction<Rank>::Ptr Psi;
	typename Wavefunction<Rank>::Ptr TempPsi;
	object OperatorCallback;

public:
	GmresWrapper() {}
	virtual ~GmresWrapper() {}

	void ApplyConfigSection(const ConfigSection &config);
	void Setup(const typename Wavefunction<Rank>::Ptr psi);
	void Solve(object callback, typename Wavefunction<Rank>::Ptr psi, typename Wavefunction<Rank>::Ptr tempPsi, bool usePypropIntegration);

	void ApplyOperator(blitz::Array<cplx, 1> &input, blitz::Array<cplx, 1> &output);

	/*
	 * Returns the _converged_ eigenvectors, a N by M matrix, where
	 * N == GetEigenvalueCount(), and M is the size of the wavefunction
	 */
	double EstimateMemoryUsage(int matrixSize, int basisSize)
	{
		return Solver.EstimateMemoryUsage(matrixSize, basisSize);
	}

	/* 
	 * Returns the nuber of matrix-vector operations performed
	 */
	int GetOperatorCount()
	{
		return Solver.GetOperatorCount();
	}

	/*
	 * Returns the number of re-orthogonalization steps performed
	 */
	int GetOrthogonalizationCount()
	{
		return Solver.GetOrthogonalizationCount();
	}

	void PrintStatistics()
	{
		Solver.PrintStatistics();
	}

	void ResetStatistics()
	{
		Solver.ResetStatistics();
	}

	double GetErrorEstimate()
	{
		return Solver.GetErrorEstimate();
	}
};

} // Namespace

#endif


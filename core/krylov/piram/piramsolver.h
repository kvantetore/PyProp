#ifndef PIRAMSOLVER_H
#define PIRAMSOLVER_H

#include <core/common.h>
#include <core/wavefunction.h>
#include <core/utility/boostpythonhack.h>

#include "piram/piram.h"

namespace krylov
{

template<int Rank>
class PiramSolver 
{
public:
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;

private:
	piram::pIRAM<cplx> Solver;

	//Temporary member variables
	Wavefunction<Rank> *Psi;
	Wavefunction<Rank> *TempPsi;
	object Callback;

	//test
	//blitz::Array<cplx, 2> M;

public:
	void ApplyConfigSection(const ConfigSection &config);
	void Setup(const Wavefunction<Rank> &psi);
	void Solve(object callback, Wavefunction<Rank> &psi, Wavefunction<Rank> &tempPsi);

	void SetupResidual(blitz::Array<cplx, 1> &residual);
	void ApplyOperator(blitz::Array<cplx, 1> &input, blitz::Array<cplx, 1> &output);

	/*
	 * Returns the converged eigenvalues
	 */
	blitz::Array<cplx, 1> GetEigenvalues()
	{
		return Solver.GetEigenvalues();
	}

	/*
	 * Returns the _converged_ eigenvectors, a N by M matrix, where
	 * N == GetEigenvalueCount(), and M is the size of the wavefunction
	 */
	blitz::Array<cplx, 1> GetEigenvector(int eigenvectorIndex)
	{
		return Solver.GetEigenvector(eigenvectorIndex);
	}

	double EstimateMemoryUsage(int matrixSize, int basisSize)
	{
		return Solver.EstimateMemoryUsage(matrixSize, basisSize);
	}

	blitz::Array<double, 1> GetErrorEstimates()
	{
		return Solver.GetErrorEstimates();
	}


	blitz::Array<double, 1> GetConvergenceEstimates()
	{
		return Solver.GetConvergenceEstimates();
	}


	/*
	 * Returns the number of converged eigenvalues
	 * <= EigenvalueCount
	 */
	int GetEigenvalueCount()
	{
		return Solver.GetConvergedEigenvalueCount();
	}
	
	/*
	 * Returns the number of restarting steps performed
	 * <= MaxIterationCount
	 */
	int GetRestartCount()
	{
		return Solver.GetRestartCount();
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

};

} // Namespace

#endif


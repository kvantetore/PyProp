#ifndef EXPOKITPROPAGATOR_H
#define EXPOKITPROPAGATOR_H

#include <core/common.h>
#include <core/wavefunction.h>
#include <core/utility/boostpythonhack.h>

#include "../krylovbase.h"

namespace krylov
{

template<int Rank>
class ArpackPropagator 
{
public:
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;

private:
    //Workspace variables
	blitz::Array<int, 1> WorkSelect;       
	blitz::Array<cplx, 1> WorkVectors;     
	blitz::Array<cplx, 1> WorkVectors2;    
	blitz::Array<cplx, 1> WorkData;        
	blitz::Array<cplx, 1> WorkResidual;    
	blitz::Array<double, 1> WorkReal;      
	blitz::Array<int, 1> IterationParameters;

	blitz::Array<cplx, 1> Eigenvalues;
	blitz::Array<cplx, 2> Eigenvectors;

public:
	//Parameters read from config file
	int BasisSize;
	double Tolerance;
	int MaxIterationCount;
	int EigenvalueCount;
	bool RandomStart;

	void ApplyConfigSection(const ConfigSection &config);
	void Setup(const Wavefunction<Rank> &psi);
	void Solve(object callback, Wavefunction<Rank> &psi, Wavefunction<Rank> &tempPsi);

	/*
	 * Returns the converged eigenvalues
	 */
	blitz::Array<cplx, 1> GetEigenvalues()
	{
		return Eigenvalues(blitz::Range(0, GetEigenvalueCount()-1));
	}

	/*
	 * Returns the _converged_ eigenvectors, a N by M matrix, where
	 * N == GetEigenvalueCount(), and M is the size of the wavefunction
	 */
	blitz::Array<cplx, 2> GetEigenvectors()
	{
		return Eigenvectors(blitz::Range(0, GetEigenvalueCount()-1), blitz::Range::all());
	}

	/*
	 * Returns the number of converged eigenvalues
	 * <= EigenvalueCount
	 */
	int GetEigenvalueCount()
	{
		return IterationParameters(5 - 1);
	}
	
	/*
	 * Returns the number of arnoldi update iterations actually performed
	 * <= MaxIterationCount
	 */
	int GetIterationCount()
	{
		return IterationParameters(3 - 1);
	}

	/* 
	 * Returns the nuber of matrix-vector operations performed
	 */
	int GetOperationCount()
	{
		return IterationParameters(9 - 1);
	}

	/*
	 * Returns the number of re-orthogonalization steps performed
	 */
	int GetOrthogonalizationCount()
	{
		return IterationParameters(11 - 1);
	}

};

} // Namespace

#endif


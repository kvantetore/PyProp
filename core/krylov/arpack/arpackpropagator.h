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

	blitz::Array<cplx, 1> Eigenvalues;
	blitz::Array<cplx, 2> Eigenvectors;

public:
	//Parameters read from config file
	int BasisSize;
	double Tolerance;
	int MaxIterationCount;
	int EigenvalueCount;

	void ApplyConfigSection(const ConfigSection &config);
	void Setup(const Wavefunction<Rank> &psi);
	void Solve(object callback, Wavefunction<Rank> &psi, Wavefunction<Rank> &tempPsi);

	blitz::Array<cplx, 1> GetEigenvalues()
	{
		return Eigenvalues;
	}

	blitz::Array<cplx, 2> GetEigenvectors()
	{
		return Eigenvectors;
	}

};

} // Namespace

#endif


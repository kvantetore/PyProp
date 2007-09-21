#ifndef EXPOKITPROPAGATOR_H
#define EXPOKITPROPAGATOR_H

#include <core/common.h>
#include <core/wavefunction.h>

#ifndef __GCCXML__
#include <boost/python.hpp>
using namespace boost::python;
#else
//UGLY hack to make pyste parse this file. 
struct object
{
//	int a;
};
#endif

namespace krylov
{

template<int Rank>
class ExpokitPropagator
{
public:
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;

private:
	blitz::Array<cplx, 1> Workspace;
	blitz::Array<int, 1> IntegerWorkspace;

public:
	Wavefunction<Rank> *Psi;
	Wavefunction<Rank> *TempPsi;
	object MultiplyCallback;
	double CurTime;
	double TimeStep;
	bool ImaginaryTime;

	int BasisSize;
	double Tolerance;
	double MatrixNorm;
	
	void ApplyConfigSection(const ConfigSection &config);
	void Setup(const Wavefunction<Rank> &psi);
	void AdvanceStep(object callback, Wavefunction<Rank> &psi, Wavefunction<Rank> &tempPsi, cplx dt, double t);

};

} // Namespace

#endif


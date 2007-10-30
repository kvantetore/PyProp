#ifndef ODEWRAPPER_H 
#define ODEWRAPPER_H

#include <core/common.h>
#include <core/wavefunction.h>

#include <core/utility/boostpythonhack.h>

namespace ODE
{

template<int Rank>
class OdeWrapper
{
public:	
	typedef blitz::Array<cplx, 1> DataArray1D;

	Wavefunction<Rank> *Psi;
	Wavefunction<Rank> *TempPsi;
	object MultiplyCallback;

private:
	double RelativeError;
	double AbsoluteError;
	int Flag;
	DataArray1D Work;
	blitz::TinyVector<int, 5> Iwork;
	double OutputTime;

public:
	void ApplyConfigSection(const ConfigSection &config);
	void Setup(const Wavefunction<Rank> &psi);
	void AdvanceStep(object callback, Wavefunction<Rank> &psi, Wavefunction<Rank> &tempPsi, cplx dt, double t);
};

}

#endif


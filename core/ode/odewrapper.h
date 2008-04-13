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

	typename Wavefunction<Rank>::Ptr Psi;
	typename Wavefunction<Rank>::Ptr TempPsi;
	object MultiplyCallback;
	bool ImTime;

private:
	double RelativeError;
	double AbsoluteError;
	int Flag;
	DataArray1D Work;
	blitz::Array<int, 1> Iwork;
	double OutputTime;
	double StartTime;

public:
	void ApplyConfigSection(const ConfigSection &config);
	void Setup(const Wavefunction<Rank> &psi);
	void AdvanceStep(object callback, typename Wavefunction<Rank>::Ptr psi, typename Wavefunction<Rank>::Ptr tempPsi, cplx dt, double t);

	double GetPropagatedTime()
	{
		return OutputTime;
	}

	void SetStartTime(double time)
	{
		this->StartTime = time;
	}
};

}

#endif


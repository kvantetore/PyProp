#ifndef EXPOKITPROPAGATOR_H
#define EXPOKITPROPAGATOR_H

#include <core/common.h>
#include <core/wavefunction.h>

#include "../krylovbase.h"
#include <core/utility/boostpythonhack.h>

namespace krylov
{

template<int Rank>
class ExpokitPropagator : public KrylovBase<Rank>
{
public:
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;

private:
	blitz::Array<cplx, 1> Workspace;
	blitz::Array<int, 1> IntegerWorkspace;

public:
	int BasisSize;
	double Tolerance;
	double MatrixNorm;
	
	void ApplyConfigSection(const ConfigSection &config);
	void Setup(const typename Wavefunction<Rank>::Ptr psi);
	void AdvanceStep(object callback, typename Wavefunction<Rank>::Ptr psi, typename Wavefunction<Rank>::Ptr tempPsi, cplx dt, double t);

};

} // Namespace

#endif


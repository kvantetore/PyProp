#ifndef KRYLOVBASE_H
#define KRYLOVBASE_H

#include <core/common.h>
#include <core/wavefunction.h>
#include <core/utility/boostpythonhack.h>

namespace krylov
{

template<int Rank>
class KrylovBase
{
public:
	Wavefunction<Rank> *Psi;
	Wavefunction<Rank> *TempPsi;
	object MultiplyCallback;
	double CurTime;
	double TimeStep;
	bool ImaginaryTime;
};

} // Namespace

#endif


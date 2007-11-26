#ifndef POTENTIALBASE_H
#define POTENTIALBASE_H

#include "../common.h"
#include "../wavefunction.h"

/* This is the base class for all potentials. Real potentials should inherit 
 * from this class in order to be forward compatible with future changes to 
 * the potential evaluation.
 * 
 * Note that we inherit from this class, but it has no virtual methods. This is
 * due to the fact that we do not need polymorphism in the normal sense, but we rather
 * use templates and inheritance to create a much faster polymorphism without vtable lookups
 *
 * Please do not use any virtual functions in your inheriting classes.
 */
template<int Rank>
class PotentialBase
{
public:
	//Updated by the PotentialEvaluator at every time step
	double CurTime;
	cplx TimeStep;

	void CurTimeUpdated() {}
	bool IsTimeDependent() { return true; }

	/* Inheriting classes should implement one of the two
	 *     double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	 *     cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	 */
};

#endif


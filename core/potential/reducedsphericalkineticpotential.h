#ifndef REDUCEDSPHERICALKINETICPOTENTIAL_H
#define REDUCEDSPHERICALKINETICPOTENTIAL_H

#include "../common.h"
#include "../wavefunction.h"
#include "dynamicpotentialevaluator.h"


/*
 * Angular Kinetic Energy potential L**2/(2 m r**2) for spherical representation
 * in spherical harmonic representation
 */
template<int Rank>
class ReducedAngularKineticEnergyPotential : public PotentialBase<Rank>
{
public:
	double Mass;
	int RadialRank;
	int LRank;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", Mass);
		config.Get("radial_rank", RadialRank);
		config.Get("l_rank", LRank);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{
		double r = pos(RadialRank);
		if (fabs(r) < 1e-10) {
			return 0;
		}
		double l = pos(LRank);
		return l*(l+1.0)/(2.0*Mass*r*r);
	}
};


#endif


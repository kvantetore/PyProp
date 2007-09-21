#ifndef SPHERICALKINETICENERGYPOTENTIAL_H
#define SPHERICALKINETICENERGYPOTENTIAL_H

#include "../common.h"
#include "../wavefunction.h"
#include "sphericaldynamicpotentialevaluator.h"
#include "dynamicpotentialevaluator.h"


/*
 * Angular Kinetic Energy potential L**2/(2 m r**2) for spherical representation
 * in spherical harmonic representation
 */
template<int Rank>
class AngularKineticEnergyPotential : public PotentialBase<Rank>
{
public:
	double Mass;
	int RadialRank;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", Mass);
		config.Get("radial_rank", RadialRank);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{
		double r = pos(RadialRank);
		if (fabs(r) < 1e-10) {
			return 0;
		}
		double l = pos(Rank-2);
		//int m = (int)pos(Rank-1);
		return l*(l+1.0)/(2.0*Mass*r*r);
	}
};

/*
 * Radial Kinetic Energy potential for the reduced wavefunction
 */
template<int Rank>
class RadialKineticEnergyPotential : public PotentialBase<Rank>
{
public:
	double Mass;
	int RadialRank;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", Mass);
		config.Get("radial_rank", RadialRank);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{
		return sqr(pos(RadialRank))/(2*Mass);
	}
};


#endif


#ifndef POLARPOTENTIAL_H
#define POLARPOTENTIAL_H

#include "dynamicpotentialevaluator.h"
#include "potentialbase.h"

template<int Rank>
class PolarKineticPotential : public PotentialBase<Rank>
{
private:
	int AngularRank;
	int RadialRank;
	double Mass;
	
public:
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("angular_rank", AngularRank);
		config.Get("radial_rank", RadialRank);
		config.Get("mass", Mass);
	}

	double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(RadialRank);
		double r2 = r * r;
		double m = pos(AngularRank);

		double V = 0;
		if (r2 > 10e-6)
		{
			V = 0.5 * (m * m) / (Mass * r2); 
		}

		return V;
	}

	bool IsTimeDependent()
	{
		return false;
	}
};


#endif


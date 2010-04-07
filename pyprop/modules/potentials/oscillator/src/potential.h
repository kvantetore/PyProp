#ifndef EXAMPLEPOTENTIALS_H
#define EXAMPLEPOTENTIALS_H

#include <core/potential/dynamicpotentialevaluator.h>


/*
 * Dynamic potential for a harmonic oscillator
 */
template<int Rank>
class HarmonicOscillatorPotential : public PotentialBase<Rank>
{
public:
	double strength;	
	
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("strength", strength);
	}
	
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double pot = 0;
		for (int i=0;i<Rank;i++)
		{
			pot += sqr(pos(i));
		}
		return 0.5 * strength * pot;
	}
};

#endif


#ifndef CARTESIANKINETICENERGYPOTENTIAL_H
#define CARTESIANKINETICENERGYPOTENTIAL_H

#include "dynamicpotentialevaluator.h"

/* 
Dynamic potential for evaluation of the kinetic energy potential for CartesianFFTEvaluator
*/
template<int Rank>
class CartesianKineticEnergyPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;
	
	void ApplyConfigSection(const ConfigSection &config)
	{
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &momentum)
	{
		double kineticPotential = 0.0;
		for (int i=0;i<Rank;i++)
		{
			kineticPotential += sqr(momentum(i));
		}
		return kineticPotential/2.0;
	}
};

#endif


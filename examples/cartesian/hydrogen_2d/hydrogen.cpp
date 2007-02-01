#include <core/wavefunction.h>
#include <core/representation/cartesianrepresentation.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class SoftColoumbPotential
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double soft;
	double charge;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("soft", soft);
		config.Get("charge", charge);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = sqr(pos(0));
		for (int i=1;i<Rank;i++)
		{
			r += sqr(pos(i));
		}
		
		return charge / sqrt(r + sqr(soft));
	}
};



#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class QDotPotential
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters

	void ApplyConfigSection(const ConfigSection &config)
	{
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		if (Rank == 1) 
		{
			double x = pos(0);
			return 0.5 * sqr(x);
		}

		if (Rank == 2)
		{
			//Coordinates
			double x = pos(0);
			double y = pos(1);

			return 0.5 * (sqr(x) + sqr(y));
		}
	}
};





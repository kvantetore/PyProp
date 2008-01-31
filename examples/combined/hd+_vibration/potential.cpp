#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class MorsePotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	double Strength;
	double Width;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("strength", Strength);
		config.Get("width", Width);
	}

	/*
	 * Called once every timestep
	 */
	void CurTimeUpdated()
	{
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(0);
		return Strength * (exp(- 2 * Width * r) - 2 * exp(- Width * r) + 1);
	}
};




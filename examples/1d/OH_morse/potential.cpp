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
	double Minimum;
	double Theta;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("strength", Strength);
		config.Get("theta", Theta);
		config.Get("minimum", Minimum);
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
		return Strength * sqr( exp(-Theta * (r - Minimum)) - 1 ) - Strength;
	}
};


template<int Rank>
class DipoleMoment : public PotentialBase<Rank>
{
public:

	//Potential parameters
	double Mu;
	double XStar;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mu", Mu);
		config.Get("x_star", XStar);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		return Mu * x * exp(-x / XStar);
	}
};


template<int Rank>
class ZhuRabitzOperator : public PotentialBase<Rank>
{
public:

	//Potential parameters
	double Gamma;
	double XMark;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("gamma", Gamma);
		config.Get("x_mark", XMark);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		return Gamma / std::sqrt(M_PI) * exp(-Gamma * Gamma * sqr(x - XMark));
	}
};

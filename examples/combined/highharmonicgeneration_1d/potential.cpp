#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>
#include <core/potential/rankonepotentialevaluator.h>

template<int Rank>
class LaserPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double PulseDuration;
	double Frequency;
	double Amplitude;

	//Calculated parameters
	double convolutionFrequency;
	double currentAmplitude;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("pulse_duration", PulseDuration);
		config.Get("frequency", Frequency);
		config.Get("amplitude", Amplitude);

		convolutionFrequency = M_PI / PulseDuration;
	}

	/*
	 * Called once every timestep
	 */
	void CurTimeUpdated()
	{
		if (CurTime > PulseDuration)
		{
			currentAmplitude = 0;
		}
		else
		{
			currentAmplitude = Amplitude;
			currentAmplitude *= sqr(sin(CurTime * convolutionFrequency));
			currentAmplitude *= cos(CurTime * Frequency);
		}
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);

		return currentAmplitude * x;
	}
};

template<int Rank>
class StarkPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double FieldStrength;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("field_strength", FieldStrength);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);

		return FieldStrength * x;
	}
};


template<int Rank>
class ModelAtomPotential : public PotentialBase<Rank>
{
public:
	//Potential parameters
	double Charge;
	double SoftParameter;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", Charge);
		config.Get("soft_param", SoftParameter);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		
		return Charge / std::sqrt(SoftParameter + x * x);
	}
};


template<int Rank>
class ModelAtomPotentialDiff : public PotentialBase<Rank>
{
public:
	//Potential parameters
	double Charge;
	double SoftParameter;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", Charge);
		config.Get("soft_param", SoftParameter);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		
		//double potential = (SoftParameter + x * x);
		//potential *= (SoftParameter + x * x);
		//potential *= (SoftParameter + x * x);
		
		return Charge * x / std::pow(SoftParameter + x*x, 3.0/2.0);
	}
};

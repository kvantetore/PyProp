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
		double z1 = std::fabs(pos(0));
		double z2 = pos(1);

		return currentAmplitude * (z1 + z2);
	}
};

template<int Rank>
class SingleActiveElectronPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Z;
	double a1;
	double a2;
	double a3;
	double a4;
	double a5;
	double a6;


	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("z", Z);
		config.Get("a1", a1);
		config.Get("a2", a2);
		config.Get("a3", a3);
		config.Get("a4", a4);
		config.Get("a5", a5);
		config.Get("a6", a6);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = std::fabs(pos(0));
		if (r < 1e-5)
		{
			return 0.0;
		}

		double V = -(Z + a1 * exp(-a2 * r) + a3 * r * exp(-a4 * r)
			+ a5 * exp(-a6 * r)) / r;

		return V;
	}

	inline double pow(double r, double n)
	{
		double r_n = 1.0;
		for (int i = 0; i < n; i++)
		{
			r_n *= r;
		}

		return r_n;
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

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("field_strength", FieldStrength);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double z1 = std::fabs(pos(0));
		double z2 = std::fabs(pos(1));

		return FieldStrength * (z1 + z2);
	}
};


template<int Rank>
class TwoElectronCorrelation1D : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Charge;
	double SoftParam;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", Charge);
		config.Get("soft_param", SoftParam);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double z1 = std::fabs(pos(0));
		double z2 = std::fabs(pos(1));
		double z12 = z1 - z2;

		return Charge / sqrt( z12 * z12 + SoftParam);
	}
};



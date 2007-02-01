#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class CoulombPotential
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double charge;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", charge);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double V;
		double r = fabs(pos(0));
		if (r < 1e-6)
		{
			V = 0.0;
		}
		else
		{
			V = charge / r;
		}
		return V;
	}
};

template<int Rank>
class PulsePotential
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Charge;
	double LaserFrequency;
	double LaserIntensity;
	double EnvelopeFrequency;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", Charge);
		config.Get("laser_frequency", LaserFrequency);
		config.Get("laser_intensity", LaserIntensity);

		double pulseDuration = 0;
		config.Get("pulse_duration", pulseDuration);
		EnvelopeFrequency = M_PI / pulseDuration;
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double t = CurTime;
		double V = 0.;
		
		//coordinates
		double r = fabs(pos(0));

		//Coulomb
		if (r < 1e-6)
		{
			V = 0.0;
		}
		else
		{
			V = Charge / r;
		}

		//Laser
        V += r * LaserIntensity * sqr(sin(EnvelopeFrequency*t)) * sin(LaserFrequency*t);

		return V;
	}
};



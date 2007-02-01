#include <core/wavefunction.h>
#include <core/potential/sphericaldynamicpotentialevaluator.h>

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
class KulanderPotential
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Charge;
	double LaserFrequency;
	double LaserIntensity;
	double LaserTurnOnTime;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", Charge);
		config.Get("laser_frequency", LaserFrequency);
		config.Get("laser_intensity", LaserIntensity);
		
		double laserTurnOnCycles = 0;
		config.Get("laser_turn_on_cycles", laserTurnOnCycles);
		LaserTurnOnTime = laserTurnOnCycles * 2 * M_PI / LaserFrequency;
		
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double t = CurTime;
		double V = 0.;
		
		//coordinates
		double r = fabs(pos(0));
		double theta = pos(1);
		//double phi = pos(2);
		double z = r * cos(theta);

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
		double turnOn = (t < LaserTurnOnTime) ? (t / LaserTurnOnTime) : 1.0;
		V += turnOn * LaserIntensity * z * sin(LaserFrequency * t);
		
		return V;
	}
};



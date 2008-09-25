#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class H2pPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Charge;
	double NuclearSeparation;
	double Softing;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", Charge);
		config.Get("softing", Softing);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		//Coordinates
		double R, r, theta;
		if (Rank == 2)
		{
			R = 2.2;
			r = std::abs(pos(0));
			theta = pos(1);
		}
		else
		{
			R = std::abs(pos(0));
			r = std::abs(pos(1));
			theta = pos(2);
		}
		
		double r2 = sqr(r) + sqr(R) / 4 + sqr(Softing);
		double z = r * cos(theta);
	
		double angDep = R * z; 
		double V1 = 1 / sqrt(r2 + angDep); 
		double V2 = 1 / sqrt(r2 - angDep); 

		return Charge * (V1 + V2);
	}
};


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
		//Coordinates
		double r, theta;
		if (Rank == 2)
		{
			r = std::abs(pos(0));
			theta = pos(1);
		}
		else
		{
			r = std::abs(pos(1));
			theta = pos(2);
		}
			
		double z = r * cos(theta);
		return currentAmplitude * z;
	}
};




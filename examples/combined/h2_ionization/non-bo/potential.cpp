#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>


template<int Rank>
class LaserPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double FieldStrength;
	double Frequency;
	double Duration;
	double PeakTime;
	double Phase;

	//Current timestep variables
	double CurrentField;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("field_strength", FieldStrength);
		config.Get("frequency", Frequency);
		config.Get("duration", Duration);
		config.Get("peak_time", PeakTime);
		config.Get("phase", Phase);
	}

	/*
	 * Called once every timestep
	 */
	void CurTimeUpdated()
	{
		double envelope = exp( - 4 * log(2) * sqr(CurTime - PeakTime) / sqr(Duration) );
		double field = FieldStrength * cos(Frequency * (CurTime - PeakTime) + Phase);

		CurrentField = field * envelope;
	}

	/*
	 * Called for every grid point, every timestep
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double z1 = pos(1);
		double z2 = pos(2);

		return CurrentField * (z1 + z2);
	}
};



template<int Rank>
class H2Potential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double NuclearSeparation;
	double NuclearSofting;
	double RepulsionSofting;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("nuclear_softing", NuclearSofting);
		config.Get("repulsion_softing", RepulsionSofting);
	}

	/*
	 * Called for every grid point.
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double R2 = pos(0)/2;

		double z1 = pos(1);
		double z2 = pos(2);
		double deltaZ = z1 - z2; 

		//Electron-Electron interaction
		double V12 = 1.0 / sqrt( sqr(deltaZ)  + sqr(RepulsionSofting));

		//Electron-Nucleus interaction
		double V1 = - 1.0 / sqrt( sqr(z1 + R2) + sqr(NuclearSofting) )
		          + - 1.0 / sqrt( sqr(z1 - R2) + sqr(NuclearSofting) );

		double V2 = - 1.0 / sqrt( sqr(z2 + R2) + sqr(NuclearSofting) )
		          + - 1.0 / sqrt( sqr(z2 - R2) + sqr(NuclearSofting) );

		//Nuclear repulsion
		double V3 = 1.0 / abs(R2);
	   
		return V1 + V2 + V12 + V3;
	}
};



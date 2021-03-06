#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>


template<int Rank>
class BornOppenheimerLaserPotential : public PotentialBase<Rank>
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
		double r = pos(0);
		double phi = pos(1);

		double z1 = r * cos(phi);
		double z2 = r * sin(phi);

		return CurrentField * (z1 + z2);
	}
};



template<int Rank>
class H2BornOppenheimerPotential : public PotentialBase<Rank>
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
		config.Get("nuclear_separation", NuclearSeparation);
		config.Get("nuclear_softing", NuclearSofting);
		config.Get("repulsion_softing", RepulsionSofting);
	}

	/*
	 * Called for every grid point.
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(0);
		double phi = pos(1);
		double R2 = NuclearSeparation / 2;

		double z1 = r * cos(phi);
		double z2 = r * sin(phi);
		double deltaZ = z1 - z2; //sqrt(2) * r * sin(phi-M_PI/4);

		//Electron-Electron interaction
		double V12 = 1.0 / sqrt( sqr(deltaZ)  + sqr(RepulsionSofting));

		//Electron-Nucleus interaction
		double V1 = - 1.0 / sqrt( sqr(z1 + R2) + sqr(NuclearSofting) )
		          + - 1.0 / sqrt( sqr(z1 - R2) + sqr(NuclearSofting) );

		double V2 = - 1.0 / sqrt( sqr(z2 + R2) + sqr(NuclearSofting) )
		          + - 1.0 / sqrt( sqr(z2 - R2) + sqr(NuclearSofting) );

		return V1 + V2 + V12;
	}
};



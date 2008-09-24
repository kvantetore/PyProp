#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class HydrogenPotential2D : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Softing;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("softing", Softing);
	}

	/*
	 * Called for every grid point.
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		if (Rank == 1)
		{
			double x = pos(0);

			return -1.0 / std::sqrt(x*x + Softing*Softing);
		}
		if (Rank == 2)
		{
			double x = pos(0);
			double y = pos(1);

			return -1.0 / std::sqrt(x*x + y*y + Softing*Softing);
		}
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
	double FieldStrength;
	double Frequency;
	double Duration;
	double PeakTime;
	double Phase;
	double PolarizationRank;

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
		config.Get("polarization_rank", PolarizationRank);
	}

	/*
	 * Called once every timestep
	 */
	void CurTimeUpdated()
	{
		double envelope = exp( - 4 * log(2) * sqr(CurTime - PeakTime) / sqr(Duration) );
		double field = FieldStrength * sin(Frequency * (CurTime - PeakTime) + Phase);

		CurrentField = field * envelope;
	}

	/*
	 * Called for every grid point, every timestep
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double z = pos(PolarizationRank);

		return CurrentField * z;
	}
};

/*
 * This is the non-dipole velocity gauge laser potential,
 *
 *     A(x,t) * p_z = A(x,t) * i * d/dz
 *
 * Since the differentiation and space dependence are in
 * seperate ranks, we can use a split-operator approach 
 * with Fourier basis. We must ensure that this potential
 * is called when the wavefunction is in a Fourier-grid
 * reprensentation.
 */
template<int Rank>
class NonDipolePotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double FieldStrength;
	double Frequency;
	double PulseStartTime;
	double Duration;
	double PeakTime;
	double Phase;
	int PolarizationRank;
	int NonDipoleRank;      // The x-rank, as in A(x,t)p_z
	double LightSpeed;

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
		config.Get("phase", Phase);
		config.Get("polarization_rank", PolarizationRank);
		config.Get("nondipole_rank", NonDipoleRank);
		config.Get("light_speed", LightSpeed);
		config.Get("pulse_start_time", PulseStartTime);
	}

	/*
	 * Called once every timestep
	 */
	void CurTimeUpdated()
	{
		//double envelope = exp( - 4 * log(2) * sqr(CurTime - PeakTime) / sqr(Duration) );
		//double field = FieldStrength * sin(Frequency * (CurTime - PeakTime) + Phase);

		//CurrentField = field * envelope;
		;
	}

	/*
	 * Called for every grid point, every timestep
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double k_z = pos(PolarizationRank);
		double x = pos(NonDipoleRank);

		double pot = 0;
		double u = CurTime - x/LightSpeed;

		//Space-time condition
		if ( (u < Duration) || (u >= PulseStartTime) )
		{	
			double envelope = sin(M_PI / Duration * u);	
			envelope *= envelope;
			envelope *= FieldStrength;
			
			pot = envelope * sin(Frequency * u + Phase);

			pot *= k_z;
		}

		return pot;
	}
};

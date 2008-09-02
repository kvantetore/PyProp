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
		double field = FieldStrength * cos(Frequency * (CurTime - PeakTime) + Phase);

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

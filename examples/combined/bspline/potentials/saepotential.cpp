#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>
#include <core/potential/rankonepotentialevaluator.h>

/*
 * Using an effective spherically symmetric potential, an outer-shell electron in 
 * in various noble gas atoms may be accurately described in situations where
 * correlations with other electrons are not important. A description of the potential
 * may be found in Tong and Lin (2005), doi:10.1088/0953-4075/38/15/001
 *
 * Parameters for Hydrogen, Helium, Neon, Argon and Argon+ are provided
 */
template<int Rank>
class SingleActiveElectronPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Effetive potential parameters
	double Z;
	double a1;
	double a2;
	double a3;
	double a4;
	double a5;
	double a6;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
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

	/*
	 * Called for every grid point at every time step. 
	 */
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
};


#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

using namespace blitz;

/*
 * Radial Kinetic Energy
 */
template<int Rank>
class KineticEnergyPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double mass;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", mass);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		return - 1. / (2. * mass);
	}

};

/* 
 * Coulomb Potential
 */
template<int Rank>
class CoupledSphericalCoulombPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	double Z;
	int RadialRank1;
	int RadialRank2;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("z", Z);
		config.Get("radial_rank1", RadialRank1);
		config.Get("radial_rank2", RadialRank2);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r1 = pos(RadialRank1);
		double r2 = pos(RadialRank2);
		return - Z/r1 - Z/r2;
	}
};



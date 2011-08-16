#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

using namespace blitz;

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


template<int Rank>
class DipoleLaserPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Charge;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", Charge);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		return -Charge * x;
	}
};


template<int Rank>
class HarmonicPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Mass;
	double Omega;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("omega", Omega);
		config.Get("mass", Mass);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		return 0.5 * Mass * Omega * x * x;
	}
};

template<int Rank>
class RegularizedCoulombPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Charge;
	double Soft;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", Charge);
		config.Get("soft_param", Soft);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		return Charge / std::sqrt(x * x + Soft * Soft);
	}
};


/*
 * Manolopoulos transmission free absorbing potential
 *
 *   See D. E. Manolopoulos, J. Chem. Phys 117, 9552, 2002.
 *
 *   energy_cutoff: wavepacket components with energy greater than this is
 *                  absorbed
 *   grid_max:      last point on the grid absorberwise
 *   delta:         accuracy parameter, determines width of absorber
 *
 */
template<int Rank>
class ManolopoulosAbsorber : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double EnergyCutoff;
	double GridMax;
	double Delta;
	double Start;

	//CAP constants
	double A;
	double B;
	double C;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("energy_cutoff", EnergyCutoff);
		config.Get("grid_max", GridMax);
		config.Get("delta", Delta);

		A = 0.112449;
		B = 0.0082735;
		C = 2.62206;

		//Calculate absorber start
		double kmin = std::sqrt(2*EnergyCutoff);
		Start = GridMax - C / (2 * Delta * kmin);
		
		cout << "Absorber starts at " << Start << std::endl;
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		cplx cap = 0;
		if (x >= Start)
		{
			double u = 2 * Delta * std::sqrt(2*EnergyCutoff) * (x - Start);
			double y = A*u - B*u*u*u + 4.0/((C-u)*(C-u)) - 4.0/((C+u)*(C+u));
			cap = -I * EnergyCutoff * y;
		}
		return cap;
	}
};


template<int Rank>
class OverlapPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	int angularRank;
	int radialRank;

	void ApplyConfigSection(const ConfigSection &config)
	{
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		return 1.0;
	}
};


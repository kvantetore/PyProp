#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>
#include <core/potential/rankonepotentialevaluator.h>

template<int Rank>
class QuantumDotPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double WellSeparation;
	double OmegaLeft;
	double OmegaRight;
	double Mass;


	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("well_separation", WellSeparation);
		config.Get("omega_left", OmegaLeft);
		config.Get("omega_right", OmegaRight);
		config.Get("mass", Mass);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);
		double V = 0;

		if (x < 0.0)
		{
			V = 0.5 * Mass * sqr(OmegaLeft) * sqr(x + WellSeparation);
		}
		else
		{
			V = 0.5 * Mass * sqr(OmegaRight) * sqr(x - WellSeparation);
		}

		return V;
	}
};

template<int Rank>
class TwoElectronCorrelation : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Charge;
	double SoftParam;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", Charge);
		config.Get("soft_param", SoftParam);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r12;

		if (Rank == 2)
		{
			double x1 = pos(0);
			double y1 = pos(1);
			double x2 = pos(2);
			double y2 = pos(3);
			r12 = sqr(x2 - x1) + sqr(y2 - y1);
		}
		else if (Rank == 4)
		{
			double x1 = pos(0);
			double x2 = pos(1);
			r12 = sqr(x2 - x1);
		}

		return Charge / sqrt( r12 + SoftParam );
	}
};



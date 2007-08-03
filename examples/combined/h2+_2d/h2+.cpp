#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class H2pPotential
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
		config.Get("nuclear_separation", NuclearSeparation);
		config.Get("softing", Softing);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		//Coordinates
		double r = fabs(pos(0));
		double theta = pos(1);

		double r2 = sqr(r) + sqr(NuclearSeparation) / 4 + sqr(Softing);
		double z = r * cos(theta);
	
		double angDep = NuclearSeparation * z; 
		double V1 = 1 / sqrt(r2 + angDep); 
		double V2 = 1 / sqrt(r2 - angDep); 

		return Charge * (V1 + V2);
	}
};


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




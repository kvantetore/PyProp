#include <core/wavefunction.h>
#include <core/potential/sphericaldynamicpotentialevaluator.h>

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
		//coordinates
		double r = fabs(pos(0));
		double theta = pos(1);
		//double phi = pos(2);
		
		double z = r * cos(theta);
		double a = sqr(r * sin(theta));
		//double x = r * sin(theta) * cos(phi);
		//double y = r * sin(theta) * sin(phi);

		//Coulomb Potential
		double V1 = Charge / sqrt(a + sqr(z + NuclearSeparation/2.0) + Softing);
		double V2 = Charge / sqrt(a + sqr(z - NuclearSeparation/2.0) + Softing);
	
		return V1 + V2;
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




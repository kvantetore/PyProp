#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class QDotPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double Omega;
	double Separation;
	double Interaction;
	double Softing;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("omega", Omega);
		config.Get("separation", Separation);
		config.Get("interaction", Interaction);
		config.Get("softing", Softing);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x1 = 0;
		double y1 = 0;
		double x2 = 0;
		double y2 = 0;

		if (Rank==4)
		{
			x1 = pos(0);
			y1 = pos(1);
			x2 = pos(2);
			y2 = pos(3);
		}
		if (Rank==2)
		{
			x1 = pos(0);
			//x2 = pos(1);
		}

		if (Rank==1)
		{
			x1 = pos(0);
		}
	
		double V1 = 0.5 * sqr(Omega) * (sqr(fabs(x1) - Separation/2) + sqr(y1));
		double V2 = 0.5 * sqr(Omega) * (sqr(fabs(x2) - Separation/2) + sqr(y2));
		double V12 = Interaction * 1 / sqrt( sqr(x1 - x2) + sqr(y1 - y2) + sqr(Softing));

		return V1 + V2 + V12;
	}
};





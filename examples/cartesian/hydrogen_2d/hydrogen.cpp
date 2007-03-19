#include <core/wavefunction.h>
#include <core/representation/cartesianrepresentation.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class SoftColoumbPotential
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double soft;
	double charge;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("soft", soft);
		config.Get("charge", charge);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = sqr(pos(0));
		for (int i=1;i<Rank;i++)
		{
			r += sqr(pos(i));
		}
		
		return charge / sqrt(r + sqr(soft));
	}
};

void testdata(int n, int count)
{
	using namespace blitz;

	Array<double, 1> r(n);
	Array<double, 2> data(n,n);
	double rmax = 100;
	double dr = rmax / n;
	r = tensor::i * dr;
	
	data = 0;
	for (int i=0; i<count; i++)
	{
		for (int x=0; x<n; x++)
		{
			for (int y=0; y<n; y++)
			{
				data(x, y) = 1 / sqrt(r(x)*r(x) + r(y) * r(y) + 0.01);
			}
		}
	}
}



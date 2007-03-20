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


#include <core/utility/blitzblas.h>

void testdata(int n, int count)
{
	using namespace blitz;

	Array<double, 1> r(n);
	Array<cplx, 2> data(n,n);
	Array<cplx, 2> data2(n,n);
	double rmax = 100;
	double dr = rmax / n;
	r = tensor::i * dr;

	data = 1;
	data2 = 0.9;
	for (int i=0; i<count; i++)
	{
		VectorElementMultiply(data, data2, data);
	}
}



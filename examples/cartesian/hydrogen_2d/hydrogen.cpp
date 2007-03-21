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

void innerprod(int n, int count)
{
	using namespace blitz;

	Array<double, 2> data(n,n);
	Array<double, 2> data2(n,n);

	data = 1;
	data2 = 0.9;
	double d = 10;
	for (int i=0; i<count; i++)
	{
		double innerprod = 0;
		for (int i=0; i<n*n; i++)
		{
			innerprod += data.data()[i] * data2.data()[i];
		}
		d /= innerprod;
	}
	cout << d << endl;
}

void innerprod2(int n, int count)
{
	using namespace blitz;

	Array<double, 2> data(n,n);
	Array<double, 2> data2(n,n);

	data = 1;
	data2 = 0.9;
	double d = 10;
	for (int i=0; i<count; i++)
	{
		d /= sum(data * data2);
	}
	cout << d << endl;
}


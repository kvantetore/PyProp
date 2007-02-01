#include <blitz/array.h>
#include <iostream>
#include <complex>

#include <core/configuration.h>

typedef double dbl;
typedef std::complex<dbl> cplx;
const cplx I(0.0, 1.0);

//#include <core/common.h>
#include <core/hltimer.h>
#include "dynamicpotentialevaluator.h"

using namespace blitz;
using namespace std;

int main(int argc, char* argv[])
{
	int size = 2000;
	
	Array<cplx, 2> array(size,size);
	array = tensor::i + tensor::j*size;

	DynamicPotentialEvaluator<CartesianKineticEnergyPotential<2>, 2> pot;

	double max = 0.0;
	for (int i=0;i<10;i++)
	{
		double gps = pot.ApplyPotential(array, 0.1, 0.0);	
		cout << "Throughput: " << gps/(1024*1024) << " Mgridpoints/second" << endl;
		if (gps > max) 
		{	
			max = gps;
		}
	}
	cout << "Max throughput: " << max/(1024*1024) << " Mgridpoints/second" << endl;
	cout << "Hei";

	return 0;
}


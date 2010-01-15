#include <math.h>
#include <stdio.h>

#include "tools.h"

// Find the weights for a "transformed" Clenshaw Curtis rule.
// written Jan. 12 2007, T.S.
// modified Jan 12 2007, torebi
//   - Modified to fit into pyprop
//   - Added scaling factors

namespace TransformedGrid 
{
		
void SetupWeights (int N, const Parameter &param, const blitz::Array<double, 1> &r, blitz::Array<double, 1> &weight)
{
	Array<double, 1> x(N);
	weight.resize(N);

	// Basic Clenshaw Curtis points and  weigths 
	x = cos(M_PI * (tensor::i +1) / (N+1));

	for (int k=0; k<N; k++)
	{
		weight(k) = 1.0;
		for (int i=1; i<(N+1)/2; i++)
		{ 
			weight(k) -= 2.0/(4*i*i-1) *cos(2.0*M_PI*i*(k+1)/(N+1));
		}
		if (N%2 == 1)
		{	
			weight(k) -= 1/(N*N+2*N) *cos((N+1)*x(N/2));
		}
		weight(k) = 2.0/(N+1) * weight(k);
	}   

	// Scaling the weights according to transformations
	if (param.Range == TransformRangeRadial)
	{
		if (param.Type == TransformTypeAlgebraic) 
		{
			weight *= param.Scaling * 2.0 / ((1+x)*(1+x));
		}
		else if (param.Type == TransformTypeTrigonometric)
		{ 
			double pi_4 = M_PI/4.0;
			weight *= param.Scaling * pi_4 / sqr(cos(pi_4*(1+x)));
		}
		else if (param.Type == TransformTypeLogarithmic) 
		{
			weight *= param.Scaling * 6.0 / (1.0 - x);
		}
		else
		{
			cout << "Invalid parameter type " << param.Type << endl;
			throw std::runtime_error("Invalid parameter type");
		}
	}	
	else if (param.Range == TransformRangeCartesian)
	{
		if (param.Type == TransformTypeAlgebraic) 
		{
			weight *= param.Scaling * 1.0 / pow(1 - x*x, 1.5);
		}
		/* Not implemented for Cartesian range
		else if (param.Type == TransformTypeTrigonometric)
		{ 
			double pi_4 = M_PI/4.0;
			weight *= 
		}
		else if (param.Type == TransformTypeLogarithmic) 
		{
			weight *= 
		}
		else
		{
			cout << "Invalid parameter type " << param.Type << endl;
			throw std::runtime_error("Invalid parameter type");
		}
		*/
	}

}
        

}; //Namespace


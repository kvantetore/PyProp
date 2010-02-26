#include "tools.h"

namespace TransformedGrid
{
/*
 Computes the diagonal scaling matrices X and Y for the
 1. and 2. derivatives of the 3 transforms
*/
void XYmat(int N, const Parameter &param, Array<double, 1> &X, Array<double, 1> &Y, Array<double, 1> &r)
{
    int k;
    Array<double, 1> xc(N);

	//resize arrays
	X.resize(N);
	Y.resize(N);
	r.resize(N);

    for (k=0; k<N; k++) 
	{
		xc(k) = cos(M_PI*(k+1)/(N+1)); 
	}

	double L = param.Scaling;
	//cout << "Transform #" << param.Type << endl;
	//cout << "Scaling    " << L << endl;

	if (param.Range == TransformRangeRadial)
	{
		if (param.Type == TransformTypeAlgebraic)
		{
		  r = L *(1 - xc) / (1 + xc);
		  X = - 2 * L / (L + r) / (L + r);
		  Y = - 2 * X / (L + r);
		}
		
		else if (param.Type == TransformTypeTrigonometric)
		{
		  r = L * tan((xc + 1.0) * M_PI/4.0);
		  X =  (4.0/M_PI) * L / (L*L + r*r);
		  Y = -(8.0/M_PI) * L * r /((L*L + r*r) * (L*L + r*r));
		} 
		
		else if (param.Type == TransformTypeLogarithmic)
		{
		  double L3 = 3 * L;
		  r = -L3*log((1-xc)/2);
		  X = 2*exp(-r/L3)/L3;
		  Y = -X/L3;
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
			r = L*xc/sqrt(1 - xc*xc);
			X = L*L /pow(L*L + r*r, 1.5);
			Y = -3*L*L*r/pow(L*L + r*r, 2.5);
		}
		/* The other param types is not implemented for a Cartesian Range
		else if (param.Type == TransformTypeTrigonometric)
		{
			r = tan(pi/4*(xc + 1));
			X = 4/pi./(1+r.*r);
			Y = -8/pi*r./(1+r.*r).^2;
		}
		else if (param.Type == TransformTypeLogarithmic)
		{
			r = -3*log((1-xc)/2);
			X = 2*exp(-r/3)/3;
			Y = -X/3;
		}
		*/
		else
		{
			cout << "Invalid parameter (type) " << param.Type << endl;
			throw std::runtime_error("Invalid parameter (type)");
		}
	}
	else
	{
			cout << "Invalid parameter (range) " << param.Type << endl;
			throw std::runtime_error("Invalid parameter (range)");
	}
    
	X = X*X;
}

};


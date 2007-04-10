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
	
	if (param.Type == Transform1)
	{
	  r = L *(1 - xc) / (1 + xc);
	  X = - 2 * L / (L + r) / (L + r);
	  Y = - 2 * X / (L + r);
    }
    
	if (param.Type == Transform2)
	{
      r = L * tan((xc + 1.0) * M_PI/4.0);
      X =  (4.0/M_PI) * L / (L*L + r*r);
      Y = -(8.0/M_PI) * L * r /((L*L + r*r) * (L*L + r*r));
	} 
    
	if (param.Type == Transform3)
	{
	  double L3 = 3 * L;
      r = -L3*log((1-xc)/2);
      X = 2*exp(-r/L3)/L3;
      Y = -X/L3;
	}
    
	X = X*X;

}

};


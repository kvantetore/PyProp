#include "tools.h"

namespace TransformedGrid
{

/* Function for setting up the Chebyshev grid and Chebyshev
   differentiation matrix */


   void cheb(int N, Array<double, 1> &x, Array<double, 2> &D)
   {
		int i,j;
		int sign;
		Array<double, 1> dsum(N+2);
		Array<double, 1> c(N+2);
		Array<double, 2> dX(N+2, N+2);
				
		//Resize output arrays
		x.resize(N+2);
		D.resize(N+2, N+2);
		
		/* Start initialize and checking input*/
		if (N==0) return;
		sign = 1;

		/* Start computing */
		for (i=0; i<=N+1; i++){
	        x(i) = cos(M_PI*i/(N+1));
			c(i) = sign; 
			sign = -sign;
		}
		c(0) = 2; c(N+1) = 2*c(N+1);
		for (i=0; i<=N+1; i++){
			for (j=0; j<=N+1; j++){
				dX(i,j) = x(i) - x(j);
				if (i==j) dX(i,j)++;
			}
		}
		for (i=0; i<=N+1; i++){
			dsum(i) = 0.0;
			for (j=0; j<=N+1; j++){
				D(i,j) = c(i)/c(j)/dX(i,j);
				dsum(i) += D(i,j);
			}
		}
		for (i=0; i<=N+1; i++) D(i,i) -= dsum(i);
	} // Done 
		
};


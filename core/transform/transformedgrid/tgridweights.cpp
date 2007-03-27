#include <math.h>
#include <stdio.h>

#include "tools.h"

// Find the weights for a "transformed" Clenshaw Curtis rule.
// written Jan. 12 2007, T.S.
// modified Jan 12 2007, torebi
//   - Modified to fit into pyprop
//   - (Added scaling factors... not yet)

namespace TransformedGrid 
{
		
void SetupWeights (int N, const Parameter &param, const blitz::Array<double, 1> &r, blitz::Array<double, 1> &weight){

     int i, k;
     double fder, pi_4;
     Array<double, 1> x(N);
     double wsum;

	 weight.resize(N);

  // Basic Clenshaw Curtis points and  weigths 
     for (k=0; k<N; k++)
        x(k) = cos(M_PI*(k+1)/(N+1));
     
     wsum = 0.0;
     for (k=0; k<N; k++){
        weight(k) = 1.0;
        for (i=1; i<(N+1)/2; i++){ 
           weight(k) -= 2.0/(4*i*i-1) *cos(2.0*M_PI*i*(k+1)/(N+1));
        }
        if (N%2 == 1) weight(k) -= 1/(N*N+2*N) *cos((N+1)*x(N/2));
        weight(k) = 2.0/(N+1) * weight(k);
     }   

  // Scaling the weights according to transformations
     if (param.Type == 1) 
        for (k=0; k<N; k++){
           fder = param.Scaling * 2.0/(1.0+x(k))/(1.0+x(k));
           weight(k) *= fder;
        }
     if (param.Type == 2){ 
        pi_4 = M_PI/4.0;
        for (k=0; k<N; k++){
           fder = param.Scaling * pi_4 / cos(pi_4*(1+x(k)))/cos(pi_4*(1+x(k)));
           weight(k) *= fder;
        }
     }
     if (param.Type == 3) 
        for (k=0; k<N; k++){
           fder = param.Scaling * 6.0/(1.0-x(k));
           weight(k) *= fder;
        }
     
}
        

}; //Namespace


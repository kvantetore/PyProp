#include "sphericalharmonics.h"

using namespace blitz;

blitz::Array<cplx, 2> CreateSphericalHarmonics(int maxL, Array<double, 2> &angularGrid)
{
	//slice theta and phi values
	Array<double, 1> thetaGrid = angularGrid(Range::all(), AngularCoordIndexTheta);
	Array<double, 1> phiGrid = angularGrid(Range::all(), AngularCoordIndexPhi);
	
	//Get associated legendre functions
	Array<cplx, 2> legendreValues = CreateAssociatedLegendre(maxL, thetaGrid);
	
	//multiply by exp( i m phi ) to get spherical harmonics
	legendreValues = legendreValues(tensor::i, tensor::j) 
			* exp( AsImaginary(GetYlmM(tensor::i) * phiGrid(tensor::j)) );
			
	return legendreValues;
}

Array<cplx, 2> CreateAssociatedLegendre(int maxL, Array<double, 1> &thetaGrid)
{
	int countL = maxL + 1;
	int countLm = sqr(countL) + 1;
	int countTheta = thetaGrid.extent(firstRank);

	Array<cplx, 2> legendreValues(countLm, countTheta);
	
	//
	// Now for the 3D ass. Jacobi polynomial (actually ass. Legendre)
	//
	for(int thetaIndex=0; thetaIndex < thetaGrid.extent(firstRank); thetaIndex++)
	{
		//double x = cos(p2grid(p));
		double x = thetaGrid(thetaIndex);
	
		for(int mIndex=0; mIndex<countL; mIndex++)
		{
			// For storing polynomial values
			double cx[countL];
	
			// Parameters dependent on mIndex
			double alpha = mIndex;
			double beta = alpha;
	
			for(int lIndex=mIndex; lIndex<countL; lIndex++)
			{
				double A,B,C;
	
				if( (lIndex-mIndex) == 0)
				{
					cx[0] = 1.0;
	
				} else if( (lIndex - mIndex) == 1){
					cx[1] = (alpha + 1.)*x;
	
				} else {
					cx[0] = 1.0;
					cx[1] = (alpha + 1.)*x;
	
					A = 2. * (lIndex-mIndex) * (lIndex+mIndex) * (2.*lIndex-2.);
					B = (2.*lIndex-2.) * (2.*lIndex-1.) * (2.*lIndex);
					C = 4.*(lIndex-1.) * (lIndex-1.) * (lIndex);
	
					cx[lIndex-mIndex] = ((B * x * cx[lIndex-mIndex-1]) - (C * cx[lIndex-mIndex-2])) / A;
				}
	
				// Calculate normalization constant
				double normDivisor =  (double)pow(2.0, 2 * mIndex + 1) * sqr(tgamma(lIndex + 1.)) ;
				double normDividend  =  (2.* lIndex + 1.) * sqr(tgamma(lIndex - mIndex + 1.));
				double norm = sqrt( normDivisor / normDividend );
	
				// Assign ass. jacobipol. value to Matrix field variable
				int lmIndex = MapYlmIndex(lIndex, mIndex);
				if(mIndex == 0)
				{
					legendreValues(lmIndex, thetaIndex) = norm * cx[lIndex-mIndex];
				}
				else
				{
					double tmp = sqrt(1. - sqr(x) );
					legendreValues(lmIndex, thetaIndex) = norm * pow(tmp,mIndex) * cx[lIndex-mIndex];
				}
			}
		}
	}
	
	return legendreValues;
}


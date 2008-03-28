// Header file for reduced spherical tools class
// Implements tools for calculating the spherical harmonics
// over the S2 sphere where phi symmetry is assumed
#ifndef REDUCED_SPHERICAL_SPHERICALTOOLS_H
#define REDUCED_SPHERICAL_SPHERICALTOOLS_H

#include "../../common.h"

//TODO: Implement support for different values of m. Right now only m=0 is implemented.

namespace ReducedSpherical
{

class ReducedSphericalTools
{
public:
	typedef shared_ptr<ReducedSphericalTools> Ptr;

private:
	blitz::Array<double, 1> ThetaGrid;		//Theta points. Note that this is only the unique theta points
	blitz::Array<double, 1> Weights;		//The weights used for integrating over omega
	
	blitz::Array<double, 2> AssocLegendrePoly; //The assosciated legendre functions evaluated in (theta) 
	   										   //This array is of the same length as weights
											   //The first rank is theta-index
											   //The second rank is l-index
	int LMax;

	void SetupQuadrature();					//Sets up the complete quadrature rules for both theta and phi
	void SetupExpansion();					//Sets up the normalized associated legendre polynomials AssocLegendrePoly
		

	void ForwardTransform_Impl(blitz::Array<cplx, 3> &input, blitz::Array<cplx, 3> &output);
	void InverseTransform_Impl(blitz::Array<cplx, 3> &input, blitz::Array<cplx, 3> &output);

public:	
	int GetLMax() { return LMax; }
	blitz::Array<double, 1> GetThetaGrid() { return ThetaGrid; }
	blitz::Array<double, 1> GetWeights() { return Weights; } 
	blitz::Array<double, 2> GetAssociatedLegendrePolynomial() { return AssocLegendrePoly; }

	void Initialize(int lmax);				//Initializes the transformation to a given lmax. If the transformation
											//has already been initalized, all allocated memory are freed and the 
											//transformation is reinitialized to the new lmax.
											
	template<int Rank> void ForwardTransform(blitz::Array<cplx, Rank> input, blitz::Array<cplx, Rank> output, int thetaRank);
											//Transforms from grid to spherical harmonic representation
											
	template<int Rank> void InverseTransform(blitz::Array<cplx, Rank> input, blitz::Array<cplx, Rank> output, int thetaRank);
											//Transforms from spherial harmonic representation back to grid
};

} //namespace

#endif


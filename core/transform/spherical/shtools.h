// shtools.h
// Header file for shtools class
// Implements tools for calculating the spherical harmonics
// over the S2 sphere
#ifndef SHTOOLS_H
#define SHTOOLS_H

#include "../../common.h"

inline int MapLmIndex(int l, int m)
{
	return l * (l + 1) + m;
}

inline double Factorial(double i)
{
	if (i <= 1) return 1;
	return i * Factorial(i-1);
}


class SphericalTransformTensorGrid
{
private:
	blitz::Array<double, 1> ThetaGrid;		//Theta points. Note that this is only the unique theta points
	blitz::Array<double, 1> PhiGrid;		//Phi points. Note that this is only the unique phi points. if 
											//you need the complete tensor product grid, use OmegaGrid instead.
	
	blitz::Array<double, 2> OmegaGrid;		//The grid points (theta, phi)
	blitz::Array<double, 1> Weights;		//The weights used for integrating over theta
	
	blitz::Array<cplx, 2> SphericalHarmonic;
	blitz::Array<cplx, 2> SphericalHarmonicDerivativeTheta;
	blitz::Array<cplx, 2> SphericalHarmonicDerivativePhi;
											//Spherical harmonic functions and its derivatives 
											//evaluated in (theta, phi), (l,m)

	blitz::Array<double, 2> AssocLegendrePoly; //The assosciated legendre functions evaluated in (theta) 
											//This array is of the same length as weights
											//The first rank is theta-index
											//The second rank is (l,m)-index
											
	int LMax;

	void SetupQuadrature();					//Sets up the complete quadrature rules for both theta and phi
	void SetupExpansion();					//Sets up the normalized associated legendre polynomials AssocLegendrePoly
		

	void ForwardTransform_Impl(blitz::Array<cplx, 3> &input, blitz::Array<cplx, 3> &output);
	void InverseTransform_Impl(blitz::Array<cplx, 3> &input, blitz::Array<cplx, 3> &output);

public:	
	static blitz::Array<double, 2> GenerateThetaQuadrature(int lmax);
											//Generates quadrature rule for theta integration. This function
											//returns a (n x 2) array, where the first column is the grid points
											//and the second column is the weights used for integration.
	static blitz::Array<double, 2> EvaluateAssociatedLegendrePolynomials(int lmax, const blitz::Array<double, 1> &theta);
											//Evaluates associated legendre polys in (l,m) and the theta points 
											//specified.

	int GetLMax() { return LMax; }
	blitz::Array<double, 1> GetThetaGrid() { return ThetaGrid; }
	blitz::Array<double, 1> GetPhiGrid() { return PhiGrid; }
	blitz::Array<double, 2> GetOmegaGrid() { return OmegaGrid; }
	blitz::Array<double, 1> GetWeights() { return Weights; } 
	blitz::Array<double, 2> GetAssociatedLegendrePolynomial() { return AssocLegendrePoly; }
	blitz::Array<cplx, 2> GetSphericalHarmonic() { return SphericalHarmonic; }
	blitz::Array<cplx, 2> GetSphericalHarmonicDerivativeTheta() { return SphericalHarmonicDerivativeTheta; }
	blitz::Array<cplx, 2> GetSphericalHarmonicDerivativePhi() { return SphericalHarmonicDerivativePhi; }

	void Initialize(int lmax);				//Initializes the transformation to a given lmax. If the transformation
											//has already been initalized, all allocated memory are freed and the 
											//transformation is reinitialized to the new lmax.
											
	template<int Rank> void ForwardTransform(blitz::Array<cplx, Rank> input, blitz::Array<cplx, Rank> output, int omegaRank);
											//Transforms from grid to spherical harmonic representation
											
	template<int Rank> void InverseTransform(blitz::Array<cplx, Rank> input, blitz::Array<cplx, Rank> output, int omegaRank);
											//Transforms from spherial harmonic representation back to grid
};

#endif


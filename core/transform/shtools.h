// shtools.h
// Header file for shtools class
// Implements tools for calculating the spherical harmonics
// over the S2 sphere
#ifndef SHTOOLS_H
#define SHTOOLS_H

#include "../common.h"

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
	blitz::Array<double, 1> Weights;		//The weights used for integrating over omega
	
	blitz::Array<double, 2> AssocLegendrePoly; //The assosciated legendre functions evaluated in (theta) 
											//This array is of the same length as weights
											//The first rank is theta-index
											//The second rank is (l,m)-index
											
	int LMax;

	void SetupQuadrature();					//Sets up the complete quadrature rules for both theta and phi
	void SetupExpansion();					//Sets up the normalized associated legendre polynomials AssocLegendrePoly
		
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

	void Initialize(int lmax);				//Initializes the transformation to a given lmax. If the transformation
											//has already been initalized, all allocated memory are freed and the 
											//transformation is reinitialized to the new lmax.
											
	void ForwardTransform(blitz::Array<cplx, 2> input, blitz::Array<cplx, 2> output);
											//Transforms from grid to spherical harmonic representation
											
	void InverseTransform(blitz::Array<cplx, 2> input, blitz::Array<cplx, 2> output);
											//Transforms from spherial harmonic representation back to grid
};


/*
 * Old routines written by eric
 */
/*
using namespace blitz;

const double PI = 3.14159265358979323846264338327950288419716939937510;

// Utility inline functions to compute the index for the different arrays
int inline MapPlmIndex(int l, int m){return l*(l+1)/2 + abs(m);}
int inline MapNlmIndex(int l, int m){return l*(l+1) + m;}
int inline MapCilmIndex(int l, int m){	return l*(l+1)+m;}
int inline MapOmegaIndex(int i, int j, int nlong){return i*nlong+j;}

// definition of SHTools class members
class SHTools{
	// private:
	public:
		bool initflag; // used to check if the variables have been initialized
		int lmax; // represents the degree of legendre polinomial
		Array<cplx,1> psi; // Values of the wavefunction
		Array<cplx,1> cilm; // Spherical Harmonics coefficients
		Array<double,1> zero; // Array of the zeros
		Array<double,1> w; // Array of the weight values
		Array<double,1> nlm; // Array of the normalization coefficients
		Array<double,2> plx; // Array of the legendre polinomials
		Array<double,1> latglq; // Array of the latitude angles in radians
		Array<double,1> longlq; // Array of the longitude angles in radians
		int norm;
		int csphase;
		int nlat; // number of angles in the latitude band
		int nlong; // number of angles in the longitude band
		double factorial(int n); //as utility function

	//public:
		//constructor and destructor
		SHTools();
		~SHTools();
		//Initialization calls
		void Initialize(const int lmax, const int norm, const int csphase);
		void Initialize(const int lmax);
		void PreGLQ();
		void PlmSchmidt(Array<double,1> p, double z);
		void PreCompute();
		void GLQGridCoord();
		void PreLegendre();
		//mutators
		void SetPsi(Array<cplx,1> psi);
		void SetCilm(Array<cplx,1> cilm);
		//extractors
		Array<double,1> GetLat();
		Array<double,1> GetLong();
		Array<double,2> GetPlx();
		Array<double,1> GetW();
		int GetNlat();
		int GetNlong();
		int GetMaxL();
		Array<cplx,1> GetPsi();
		Array<cplx,1> GetCilm();
		//MAIN FUNCTIONS
		//compute the wave function psi
		void FComputePsi(Array<cplx,2> &, Array<cplx,2> &);
		Array<cplx,1> FComputePsi();
		Array<cplx,1> ComputePsi();
		//compute the coeficients
		void FComputeCilm(Array<cplx,2> &, Array<cplx,2> &);
		Array<cplx,1> FComputeCilm();
		Array<cplx,1> ComputeCilm();
};
*/

#endif


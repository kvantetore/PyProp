// shtools.cpp file

#include <iostream>
#include "shtools.h"
#include "fftw3.h"


void SphericalTransformTensorGrid::ForwardTransform(blitz::Array<cplx, 2> input, blitz::Array<cplx, 2> output)
{
	int thetaCount = ThetaGrid.extent(0);
	int phiCount = PhiGrid.extent(0);		
	int radialCount = input.extent(0);		//The first rank is the radial index
	int omegaCount = input.extent(1);		//The second rank is the omega index

	if (omegaCount != phiCount * thetaCount)
	{
		throw std::runtime_error("Something is wrong with phicount, thetacount or omegacount");
	}

	//Reinterpret the 2D INPUT data as a 3D array, where the compressed rank omega
	//is expanded into theta and phi.
	//rank0 = radius
	//rank1 = theta
	//rank2 = phi
	blitz::TinyVector<int, 3> shape(radialCount, thetaCount, phiCount);
	blitz::TinyVector<int, 3> stride(omegaCount, phiCount, 1);
	blitz::Array<cplx, 3> expandedInput(input.data(), shape, stride, blitz::neverDeleteData);
	
	// Creates a plan to perform FFT transform. 
	fftw_complex *data = (fftw_complex *) expandedInput.data();
	int transformCount = radialCount * thetaCount;
	fftw_plan plan = fftw_plan_many_dft(
						1, 				//Rank of the transform
						&phiCount,		//Length of the vector 
						transformCount,	//Number of vectors to transform
						data, 			//Input data
						NULL,		 	//no padding (input)
						1,				//Stride of vector (input)
						phiCount, 		//Distance between vectors (input)
						data, 			//Output data
						NULL,			//no padding (output)
						1,				//Stride of vector (output) 
						phiCount,		//Distance between vectors (output)
						FFTW_FORWARD,   //fft-direction 
						FFTW_MEASURE   //using estimate
						);
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	output = 0.0;
	for (int r=0; r<radialCount; r++)
	{
		for (int i=0; i<thetaCount; i++)
		{
			int lmIndex = 0;
			for (int l=0; l <=LMax ; l++)
			{
				for (int m=-l; m<0; m++)
				{
					output(r, lmIndex) += expandedInput(r, i, phiCount + m) * AssocLegendrePoly(i, lmIndex) *  Weights(i); 
					lmIndex++;
				}
				for (int m=0; m<=l; m++)
				{
					output(r, lmIndex) += expandedInput(r, i, m) * AssocLegendrePoly(i, lmIndex) *  Weights(i); 
					lmIndex++;
				}
			}
		}
	}
}

void SphericalTransformTensorGrid::InverseTransform(blitz::Array<cplx, 2> input, blitz::Array<cplx, 2> output)
{
	int thetaCount = ThetaGrid.extent(0);
	int phiCount = PhiGrid.extent(0);		
	int radialCount = input.extent(0);		//The first rank is the radial index
	int omegaCount = output.extent(1);		//The second rank is the omega index

	if (omegaCount != phiCount * thetaCount)
	{
		throw std::runtime_error("Something is wrong with phicount, thetacount or omegacount");
	}

	//Reinterpret the 2D OUTPUT data as a 3D array, where the compressed rank omega
	//is expanded into theta and phi.
	//rank0 = radius
	//rank1 = theta
	//rank2 = phi
	blitz::TinyVector<int, 3> shape(radialCount, thetaCount, phiCount);
	blitz::TinyVector<int, 3> stride(omegaCount, phiCount, 1);
	blitz::Array<cplx, 3> expandedOutput(output.data(), shape, stride, blitz::neverDeleteData);
	
	//initialize output to zero
	expandedOutput = 0;

	for (int rd=0; rd<radialCount; rd++)
	{
		for(int i=0; i<thetaCount; i++)
		{
			//iterate over (l,m) in principal order
			int lmIndex = 0;
			for (int l=0; l <=LMax; l++) 
			{
				//negative m
				for (int m=-l; m < 0; m++)
				{
					expandedOutput(rd, i, phiCount + m) += input(rd,lmIndex) * AssocLegendrePoly(i,lmIndex);
					lmIndex++;
				}
				
				//positive m
				for (int m=0; m<=l; m++)
				{
					expandedOutput(rd, i, m) += input(rd,lmIndex) * AssocLegendrePoly(i,lmIndex);
					lmIndex++;
				}
			}
		}
	}

	// Creates a plan to perform FFT transform. 
	fftw_complex *data = (fftw_complex *) expandedOutput.data();
	int transformCount = radialCount * thetaCount;
	fftw_plan plan = fftw_plan_many_dft(
						1, 				//Rank of the transform
						&phiCount,		//Length of the vector 
						transformCount,	//Number of vectors to transform
						data, 			//Input data
						NULL,		 	//no padding (input)
						1,				//Stride of vector (input)
						phiCount, 		//Distance between vectors (input)
						data, 			//Output data
						NULL,			//no padding (output)
						1,				//Stride of vector (output) 
						phiCount,		//Distance between vectors (output)
						FFTW_BACKWARD,  //fft-direction 
						FFTW_MEASURE   //using estimate
						);
	

	fftw_execute(plan);
	fftw_destroy_plan(plan);
}

/*
 * Initializes the transform to a given lmax
 */
void SphericalTransformTensorGrid::Initialize(int lmax)
{
	LMax = lmax;
	SetupQuadrature();
	SetupExpansion();
}

/*
 * Sets up the quadrature rules for integrating polynomials on the sphere
 * This function allocates and initializes ThetaGrid, PhiGrid, OmegaGrid 
 * and Weights according to LMax
 */
void SphericalTransformTensorGrid::SetupQuadrature()
{
	//In order to do integration correct up to l=lmax, we use a grid
	//of (2*lmax +1, 2*lmax +1) for our grid.
	int thetaCount = 2 * LMax + 1;
	int phiCount = 2 * LMax + 1; 
	int gridSize = thetaCount * phiCount;

	//Resize grid arrays
	ThetaGrid.resize(thetaCount);
	PhiGrid.resize(phiCount);
	OmegaGrid.resize(gridSize,2);
	Weights.resize(thetaCount);
	
	//Find theta gridpoints and weights.
	blitz::Array<double, 2> thetaQuadrature = GenerateThetaQuadrature(LMax);

	//Set up theta and phi grid
	ThetaGrid = thetaQuadrature(blitz::Range::all(), 0);
	double deltaPhi = 2.0 * M_PI / (double)phiCount;
	PhiGrid = deltaPhi * blitz::tensor::i;

	//Set up weights
	Weights = thetaQuadrature(blitz::Range::all(), 1) * deltaPhi;

	//Set up full omega grid 
	int omegaIndex = 0;
	for (int thetaIndex = 0; thetaIndex < thetaCount; thetaIndex++)
	{
		for (int phiIndex = 0; phiIndex < phiCount; phiIndex++)
		{
			OmegaGrid(omegaIndex, 0) = ThetaGrid(thetaIndex);
			OmegaGrid(omegaIndex, 1) = PhiGrid(phiIndex);
			omegaIndex++;
		}
	}
}

/* Sets up the expansion from spherical harmonics to grid basis.
 * This involves evaluating the spherical harmonic in the grid points set up
 * by SetupQuadrature() (which must be called before this function)
 * This function initializes the arrays AssocLegendePoly
 */
void SphericalTransformTensorGrid::SetupExpansion()
{
	std::cout << "Setting up expansion" << std::endl;

	blitz::Array<double, 2> legendre = EvaluateAssociatedLegendrePolynomials(LMax, ThetaGrid);
	AssocLegendrePoly.reference(legendre);

	//Normalize associated legendre such that the int(legendre_lm * legendre_lm', omega) = delta_lm_lm'
	int lmIndex = 0;
	for (int l=0; l<=LMax; l++)
	{
		for (int m=-l; m<=l; m++)
		{
			//Calculate norm(l, m)
		    double norm;
			double sign = 1 ;// (m > 0) ? pow(-1., m) : 1;
			norm = sign * sqrt(((2.0*l+1)/(4*M_PI))*(Factorial(l-abs(m))/Factorial(l+abs(m))));
			norm = sign * sqrt(((2.0*l+1)/(4*M_PI))*(Factorial(l-abs(m))/Factorial(l+abs(m))));

			//Scale assoc legendre
			AssocLegendrePoly(blitz::Range::all(), lmIndex) *= norm;

			//we are iterating (l,m) in principal order.
			lmIndex++;	
		}
	}
}


/* This function generates a quadrature for integrating polynomials up to order lmax
 * with legendre weight function. The grid points used for this integration
 * is the zeros of the (lmax+1)th order polynomial.
 *
 * Tore, Nov 10 2006: I'm not completely sure how this function works, or if it does
 * exaclty what it states above. The function is reimplemented from the fortran library
 * SHTools by Eric and it seems like it is working, so please don't touch if you don't really
 * have to...
 */
blitz::Array<double, 2> SphericalTransformTensorGrid::GenerateThetaQuadrature(int lmax)
{
	int iter;
	double p1,p2,p3,pp,z,z1;
	int m= lmax+1;

	blitz::Array<double, 2> quadrature(2 * lmax + 1, 2);
	quadrature = 0;
		
	for(int i=0; i<m; i++){
		iter = 0;
		// Approximation for the ith root
		z = cos(M_PI*(i+0.75)/(2*lmax +1.5));
		
		// Find the true value using newtons method
		while(true){
			++iter;
			p1 = 1.0;
			p2 = 0.0;
			// determine the legendre polynomial evaluated at z(p1)
			// using recurrence relationships
			for (int j=1; j<=2*lmax+1;j++){
				p3=p2;
				p2=p1;
				p1=((double)(2*j-1)*z*p2-(double)(j-1)*p3)/(double)j;
			}
			// This is the derivative of the legendre polynomial
			// using recurrence relationships
			// Newton Method
			pp = ((double)(2*lmax+1)*(z*p1-p2))/(z*z-1.0);

			z1=z;
			z = z1 - p1/pp;

			
			if (fabs(z-z1) <= 1e-14/*0.00000000000000000001*/)
				break;

			// prevent infinite loop
			if (iter > 100000)
				break;
			
		}

		if (iter > 100000)
		{
			std::cout << "Could not accurately compute zeros of the " << i << "th order legendre polynomial" << std::endl;
			std::cout << "Relative error: " << fabs(z - z1) / fabs(z) << std::endl;
			throw std::runtime_error("Error computing zeros and weights of the legendre polynomial");
		}

		//We're using twice as many points as for a real gaussian quadrature
		int i2 = 2*lmax - i;
		
		//Grid point:
		quadrature(i , 0) = acos( z); //asin( z) - M_PI / 2.0; //HMM: asin why not acos?
		quadrature(i2, 0) = acos(-z); //asin(-z) - M_PI / 2.0;

		//Weights
		quadrature(i , 1) = 2.0 / ( (1.0 - sqr(z)) * sqr(pp));
		quadrature(i2, 1) = quadrature(i, 1);
	}

	return quadrature;
}

/*
 * Evaluates the Associated Legendre Polynomial Plm(cos(theta)) for all (l,m) such that
 * l <= lmax in the theta points specified.
 *
 * Tore Nov 10 2006: This function is reimplemented from the fortran library SHTools by Eric.
 * I do not completely understand what is does and why it works, but it looks like it's computing
 * the assoc legendre poly. It is sparesely commented and using cryptic variable names, so i reccomend 
 * to avoid doing anything here unless you really have to.
 */ 
blitz::Array<double, 2> SphericalTransformTensorGrid::EvaluateAssociatedLegendrePolynomials(int lmax, const blitz::Array<double, 1> &theta)
{
	double x, somx2, fact, pmm, pmmp1,pll;
	blitz::Array<double,1> pmml(lmax+1);
	
	int thetaCount = theta.extent(0);
	int lmCount = sqr(lmax + 1);

	blitz::Array<double, 2> legendre(thetaCount, lmCount);
	

	for(int i = 0; i < thetaCount; i++){
		x = cos(theta(i)); // x=cos(theta')=cos(latglq)
		somx2 = sqrt((1.0-x)*(1.0+x));

		// First: calculate the assoc legendre for abs(m) = l
		for (int l=0; l <=lmax; l++){
			fact = 1.0;
			pmml(l) = 1.0;
			for (int k = 0; k < l; k++){
				pmml(l) *= (-fact*somx2);
				fact += 2.0;
			}
			legendre(i, MapLmIndex(l, l)) = pmml(l);
			legendre(i, MapLmIndex(l,-l)) = pmml(l);
		}

		// Second: calculate assoc legendre for abs(m) < l
		for( int l=1; l <=lmax ; l++){
			for(int m=0; m<l; m++){
				pmm = pmml(m);
				pmmp1 = x*(2.0*m+1.0)*pmm;
				if (l != (m+1)){
					for (int ll = m+2; ll< l+1; ll++){
						pll = (x*(2.0*ll-1.0)*pmmp1 - (ll+m-1.0)*pmm)/(ll-m);
						pmm = pmmp1;
						pmmp1 = pll;
					}
				}
				//Except from a normalization constant, the associated legendre
				//poly are symmetric for m -> -m
				legendre(i, MapLmIndex(l, m)) = pmmp1;
				legendre(i, MapLmIndex(l,-m)) = pmmp1;
			}
		}
	}

	return legendre;
}



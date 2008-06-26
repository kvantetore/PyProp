#include "reducedsphericaltools.h"
#include "../shtools.h"
#include "../../utility/blitztricks.h"
#include "../../utility/blitzblas.h"
#include "../../utility/orthogonalpolynomial.h"

namespace ReducedSpherical
{

int ReducedSphericalTools::GetAlgorithm(int preCout, int postCount)
{
	int algo = Algorithm;
	if (algo == -1)
	{
		if (postCount == 1)
		{
			algo = 3;
		}
		else
		{
			algo = 0;
		}
	}
	return algo;
}

template<int Rank> void ReducedSphericalTools::ForwardTransform(blitz::Array<cplx, Rank> input, blitz::Array<cplx, Rank> output, int thetaRank)
{
	//Project input and output into 3D arrays with thetaRank in the middle. if omega is highest 
	//or lowest rank, the highest or lowest rank will be of size 1 
	blitz::Array<cplx, 3> input3d = MapToRank3(input, thetaRank, 1);
	blitz::Array<cplx, 3> output3d = MapToRank3(output, thetaRank, 1);

	ForwardTransform_Impl(input3d, output3d);
}
	
template<int Rank> void ReducedSphericalTools::InverseTransform(blitz::Array<cplx, Rank> input, blitz::Array<cplx, Rank> output, int thetaRank)
{
	//Project input and output into 3D arrays with thetaRank in the middle. if omega is highest 
	//or lowest rank, the highest or lowest rank will be of size 1 
	blitz::Array<cplx, 3> input3d = MapToRank3(input, thetaRank, 1);
	blitz::Array<cplx, 3> output3d = MapToRank3(output, thetaRank, 1);

	InverseTransform_Impl(input3d, output3d);
}

void ReducedSphericalTools::ForwardTransform_Impl(blitz::Array<cplx, 3> &input, blitz::Array<cplx, 3> &output)
{
	int preCount = input.extent(0);
	int postCount = input.extent(2);
	int thetaCount = ThetaGrid.extent(0);

	int algo = GetAlgorithm(preCount, postCount);

	if (algo == 0)
	{
		output = 0;
		for (int i=0; i<preCount; i++)
		{
			for (int lIndex=0; lIndex<thetaCount; lIndex++)
			{
				for (int thetaIndex=0; thetaIndex<thetaCount; thetaIndex++)
				{
					double legendre = AssocLegendrePolyTransposed(lIndex, thetaIndex) * Weights(thetaIndex);
		
					for (int j=0; j<postCount; j++)
					{
						output(i, lIndex, j) +=  legendre * input(i, thetaIndex, j);
					}
				}
			}
		}
	}
	else if (algo == 1)
	{
		blitz::Array<cplx, 2> input2d = MapToRank2(input, 1);
		blitz::Array<cplx, 2> output2d = MapToRank2(output, 1);

		for (int i=0; i<preCount; i++)
		{
			for (int lIndex=0; lIndex<thetaCount; lIndex++)
			{
				output2d(i, lIndex) = 0;
				for (int thetaIndex=0; thetaIndex<thetaCount; thetaIndex++)
				{
					double legendre = AssocLegendrePolyTransposed(lIndex, thetaIndex) * Weights(thetaIndex);
					output2d(i, lIndex) +=  legendre * input(i, thetaIndex);
				}
			}
		}
	}
	else if (algo == 2)
	{
		blitz::Array<cplx, 2> input2d = MapToRank2(input, 1);
		blitz::Array<cplx, 2> output2d = MapToRank2(output, 1);
	
		cplx* __restrict__ outputPtr = &output2d(0, 0);
		
		cplx* __restrict__ inputPtr;
		cplx* __restrict__ inputStartPtr = &input2d(0, 0);
		
		double* __restrict__ legendrePtr;
		double* __restrict__ legendreStartPtr = &AssocLegendrePolyTransposed(0, 0);

		double* __restrict__ weightsPtr;
		double* __restrict__ weightsStartPtr = &Weights(0,0);

		for (int i=0; i<preCount; i++)
		{
			legendrePtr = legendreStartPtr;
			for (int lIndex=0; lIndex<thetaCount; lIndex++)
			{
				inputPtr = inputStartPtr;
				weightsPtr = weightsStartPtr;

				*outputPtr = 0;
				for (int thetaIndex=0; thetaIndex<thetaCount; thetaIndex++)
				{
					double legendre = (*legendrePtr++) * (*weightsPtr++);
					*outputPtr +=  legendre * (*inputPtr++);
				}
				outputPtr++;
			}
			inputStartPtr+= input2d.stride(0);
		}
	}
	else if (algo == 3)
	{
		blitz::Array<cplx, 2> input2d = MapToRank2(input, 1);
		blitz::Array<cplx, 2> output2d = MapToRank2(output, 1);		

		for (int i=0; i<preCount; i++)
		{
			blitz::Array<cplx, 1> inputSlice = input2d(i, blitz::Range::all());
			blitz::Array<cplx, 1> outputSlice = output2d(i, blitz::Range::all());

			MatrixVectorMultiply(ForwardMatrix, inputSlice, outputSlice);
		}
	}
	else if (algo == 4)
	{
		for (int i=0; i<preCount; i++)
		{
			blitz::Array<cplx, 2> inputSlice = input(i, blitz::Range::all(), blitz::Range::all());
			blitz::Array<cplx, 2> outputSlice = output(i, blitz::Range::all(), blitz::Range::all());

			MatrixMatrixMultiply(ForwardMatrix, inputSlice, outputSlice);
		}
	}
}

void ReducedSphericalTools::InverseTransform_Impl(blitz::Array<cplx, 3> &input, blitz::Array<cplx, 3> &output)
{
	int thetaCount = ThetaGrid.extent(0);
	int preCount = input.extent(0);
	int postCount = input.extent(2);

	int algo = GetAlgorithm(preCount, postCount);

	if (algo == 0)
	{
		output = 0;
		for (int i=0; i<preCount; i++)
		{
			for (int thetaIndex=0; thetaIndex<thetaCount; thetaIndex++)
			{
				for (int lIndex=0; lIndex<thetaCount; lIndex++)
				{
					double legendre = AssocLegendrePoly(thetaIndex, lIndex);
		
					for (int j=0; j<postCount; j++)
					{
						output(i, thetaIndex, j) +=  legendre * input(i, lIndex, j);
					}
				}
			}
		}
	}

	else if (algo == 1)
	{
		blitz::Array<cplx, 2> input2d = MapToRank2(input, 1);
		blitz::Array<cplx, 2> output2d = MapToRank2(output, 1);
	
		for (int i=0; i<preCount; i++)
		{
			for (int thetaIndex=0; thetaIndex<thetaCount; thetaIndex++)
			{
				output2d(i, thetaIndex) = 0;
				for (int lIndex=0; lIndex<thetaCount; lIndex++)
				{
					double legendre = AssocLegendrePoly(thetaIndex, lIndex);
					output2d(i, thetaIndex) +=  legendre * input2d(i, lIndex);
				}
			}
		}
	}

	else if (algo == 2)
	{
		blitz::Array<cplx, 2> input2d = MapToRank2(input, 1);
		blitz::Array<cplx, 2> output2d = MapToRank2(output, 1);
	
		cplx* __restrict__ outputPtr;
		cplx* __restrict__ outputStartPtr = &output2d(0, 0);
		
		cplx* __restrict__ inputPtr;
		cplx* __restrict__ inputStartPtr = &input2d(0, 0);
		
		double* __restrict__ legendrePtr;
		double* __restrict__ legendreStartPtr = &AssocLegendrePoly(0, 0);

		outputPtr = outputStartPtr;
		for (int i=0; i<preCount; i++)
		{
			legendrePtr = legendreStartPtr;
			for (int thetaIndex=0; thetaIndex<thetaCount; thetaIndex++)
			{
				*outputPtr = 0;
				inputPtr = inputStartPtr;
				for (int lIndex=0; lIndex<thetaCount; lIndex++)
				{
					double legendre = (*legendrePtr++);
					*outputPtr +=  legendre * (*inputPtr++);
				}
				outputPtr++;
			}
			inputStartPtr += input2d.stride(0);
		}
	}

	else if (algo == 3)
	{
		blitz::Array<cplx, 2> input2d = MapToRank2(input, 1);
		blitz::Array<cplx, 2> output2d = MapToRank2(output, 1);		

		for (int i=0; i<preCount; i++)
		{
			blitz::Array<cplx, 1> inputSlice = input2d(i, blitz::Range::all());
			blitz::Array<cplx, 1> outputSlice = output2d(i, blitz::Range::all());

			MatrixVectorMultiply(InverseMatrix, inputSlice, outputSlice);
		}
	}

	else if (algo == 4)
	{
		for (int i=0; i<preCount; i++)
		{
			blitz::Array<cplx, 2> inputSlice = input(i, blitz::Range::all(), blitz::Range::all());
			blitz::Array<cplx, 2> outputSlice = output(i, blitz::Range::all(), blitz::Range::all());

			MatrixMatrixMultiply(InverseMatrix, inputSlice, outputSlice);
		}
	}

}

/*
 * Initializes the transform to a given lmax
 */
void ReducedSphericalTools::Initialize(int lmax)
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
void ReducedSphericalTools::SetupQuadrature()
{
	//In order to do integration correct up to l=lmax, we use a grid
	//of lmax+1 collocation points.
	int thetaCount = LMax + 1;

	//Resize grid arrays
	ThetaGrid.resize(thetaCount);
	Weights.resize(thetaCount);
	
	//Find theta gridpoints and weights.
	blitz::Array<double, 2> thetaQuadrature = OrthogonalPolynomial::CalculateQuadrature(LMax+1, OrthogonalPolynomial::Legendre1);

	//Set up theta grid
	ThetaGrid = acos( thetaQuadrature(blitz::Range::all(), 0) );

	//Set up weights
	Weights = thetaQuadrature(blitz::Range::all(), 1);
}

/* Sets up the expansion from spherical harmonics to grid basis.
 * This involves evaluating the spherical harmonic in the grid points set up
 * by SetupQuadrature() (which must be called before this function)
 * This function initializes the arrays AssocLegendePoly
 */
void ReducedSphericalTools::SetupExpansion()
{
	//TODO: Add support for different fixed m's
	std::cout << "Setting up expansion" << std::endl;

	int m = 0; //fixed m

	int thetaCount = ThetaGrid.extent(0);
	//Calculate all assoc legendre poly, even though we only need some of them. and why? because it works!
	blitz::Array<double, 2> legendre = SphericalTransformTensorGrid::EvaluateAssociatedLegendrePolynomials(LMax, ThetaGrid);
	AssocLegendrePoly.resize(thetaCount, thetaCount); //use as many theta points as sph harm basis funcs.
	AssocLegendrePolyTransposed.resize(thetaCount, thetaCount); //use as many theta points as sph harm basis funcs.
	AssocLegendrePolyDerivative.resize(thetaCount, thetaCount);

	InverseMatrix.resize(thetaCount, thetaCount);
	ForwardMatrix.resize(thetaCount, thetaCount);

	blitz::Range all = blitz::Range::all();

	//Normalize associated legendre such that the int(legendre_lm * legendre_lm', omega) = delta_lm_lm'
	for (int l=0; l<=LMax; l++)
	{
		//Calculate norm(l, m)
	    double norm;
		double sign = 1 ;// (m > 0) ? pow(-1., m) : 1;
		norm = sign * sqrt(((2.0*l+1)/(4*M_PI))*(Factorial(l-abs(m))/Factorial(l+abs(m))));
		norm = sign * sqrt(((2.0*l+1)/(4*M_PI))*(Factorial(l-abs(m))/Factorial(l+abs(m))));
		norm = norm * sqrt(2*M_PI);

		//Scale assoc legendre
		AssocLegendrePoly(blitz::Range::all(), l) = legendre(blitz::Range::all(), MapLmIndex(l, m)) *  norm;
		AssocLegendrePolyTransposed(l, blitz::Range::all()) = legendre(blitz::Range::all(), MapLmIndex(l, m)) *  norm;


		AssocLegendrePolyDerivative(blitz::Range::all(), l) = 0;
		if (l > 0)
		{
			double normalization = sqrt( (2*l+1.)/(2*l-1.) * (l-std::abs(m)) / (double)(l+std::abs(m)));
			AssocLegendrePolyDerivative(all, l) = l * cos(ThetaGrid) * AssocLegendrePoly(all, l) - (l+m) * normalization * AssocLegendrePoly(all,l-1);
			//AssocLegendrePolyDerivative(all, l) /= sin(ThetaGrid);
		}


		for (int j=0; j<thetaCount; j++)
		{
			ForwardMatrix(l, j) = legendre(j, MapLmIndex(l, m)) * norm * Weights(j);
			InverseMatrix(j, l) = legendre(j, MapLmIndex(l, m)) * norm;
		}
	}
}



template void ReducedSphericalTools::ForwardTransform(blitz::Array<cplx, 1> input, blitz::Array<cplx, 1> output, int thetaRank);
template void ReducedSphericalTools::ForwardTransform(blitz::Array<cplx, 2> input, blitz::Array<cplx, 2> output, int thetaRank);
template void ReducedSphericalTools::ForwardTransform(blitz::Array<cplx, 3> input, blitz::Array<cplx, 3> output, int thetaRank);

template void ReducedSphericalTools::InverseTransform(blitz::Array<cplx, 1> input, blitz::Array<cplx, 1> output, int thetaRank);
template void ReducedSphericalTools::InverseTransform(blitz::Array<cplx, 2> input, blitz::Array<cplx, 2> output, int thetaRank);
template void ReducedSphericalTools::InverseTransform(blitz::Array<cplx, 3> input, blitz::Array<cplx, 3> output, int thetaRank);

} //namespace


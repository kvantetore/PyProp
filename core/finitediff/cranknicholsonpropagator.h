#ifndef CRANKNICHOLSONPROPAGATOR_H
#define CRANKNICHOLSONPROPAGATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/combinedrepresentation.h"
#include "../representation/cartesianrepresentation.h"
#include "../utility/blitzlapack.h"
#include "../utility/blitzblas.h"
#include "../utility/blitztricks.h"

/*
 * Propagates the 1-rank kinetic energy and an optional time-independent potential 
 * using the crank nicholson scheme.
 *
 * Any grid may be used, as long as it is ordered, that is x_i < x_i+1 for all i
 *
 * Finite differences are used to approximate the 1-rank kinetic energy operator 1/(2 m) d2/dx2
 *
 */

/* 
 * These are defined in combinedrepresentation_generated.cpp
 */
//Rank1 Hermitian Banded
void TensorPotentialMultiply_Rank1_Band(int rank, blitz::Array<cplx, 1> potential, double scaling, blitz::Array<cplx, 1> &source, blitz::Array<cplx, 1> &dest);
void TensorPotentialMultiply_Rank1_Band(int rank, blitz::Array<cplx, 2> potential, double scaling, blitz::Array<cplx, 2> &source, blitz::Array<cplx, 2> &dest);
void TensorPotentialMultiply_Rank1_Band(int rank, blitz::Array<cplx, 3> potential, double scaling, blitz::Array<cplx, 3> &source, blitz::Array<cplx, 3> &dest);
void TensorPotentialMultiply_Rank1_Band(int rank, blitz::Array<cplx, 4> potential, double scaling, blitz::Array<cplx, 4> &source, blitz::Array<cplx, 4> &dest);

//Rank1 Non-Hermitian Banded
void TensorPotentialMultiply_Rank1_BandNH(int rank, blitz::Array<cplx, 1> potential, double scaling, blitz::Array<cplx, 1> &source, blitz::Array<cplx, 1> &dest);
void TensorPotentialMultiply_Rank1_BandNH(int rank, blitz::Array<cplx, 2> potential, double scaling, blitz::Array<cplx, 2> &source, blitz::Array<cplx, 2> &dest);
void TensorPotentialMultiply_Rank1_BandNH(int rank, blitz::Array<cplx, 3> potential, double scaling, blitz::Array<cplx, 3> &source, blitz::Array<cplx, 3> &dest);
void TensorPotentialMultiply_Rank1_BandNH(int rank, blitz::Array<cplx, 4> potential, double scaling, blitz::Array<cplx, 4> &source, blitz::Array<cplx, 4> &dest);


//Rank1 Non-Hermitian Banded
void TensorPotentialMultiply_Rank1_Dense(int rank, blitz::Array<cplx, 1> potential, double scaling, blitz::Array<cplx, 1> &source, blitz::Array<cplx, 1> &dest);
void TensorPotentialMultiply_Rank1_Dense(int rank, blitz::Array<cplx, 2> potential, double scaling, blitz::Array<cplx, 2> &source, blitz::Array<cplx, 2> &dest);
void TensorPotentialMultiply_Rank1_Dense(int rank, blitz::Array<cplx, 3> potential, double scaling, blitz::Array<cplx, 3> &source, blitz::Array<cplx, 3> &dest);
void TensorPotentialMultiply_Rank1_Dense(int rank, blitz::Array<cplx, 4> potential, double scaling, blitz::Array<cplx, 4> &source, blitz::Array<cplx, 4> &dest);



template<int Rank>
class CrankNicholsonPropagator
{
public:

	/* 
	 * Multiplies the kinetic energy operator on sourcePsi and 
	 * adds the result to destPsi.
	 *
	 * Note that it does not multiply the potential
	 */
	void MultiplyKineticEnergyOperator(Wavefunction<Rank> &sourcePsi, Wavefunction<Rank> &destPsi)
	{
		//Reshape the kinetic energy matrix into a N-d array suitable for TensorPotentialMultiply
		blitz::TinyVector<int, Rank> shape = 1;
		shape(TransformRank) = LaplacianBlasBanded.size();
		blitz::TinyVector<int, Rank> stride = 1;
		for (int j=0; j<TransformRank; j++)
		{
			stride(j) = LaplacianBlasBanded.size();
		}

		//Use TensorPotential mechanism to perform multiple matrix-vector multiplications
		blitz::Array<cplx, Rank> src = sourcePsi.GetData();
		blitz::Array<cplx, Rank> dst = destPsi.GetData();
		
		blitz::Array<cplx, Rank> kineticEnergyTensor(LaplacianBlasBanded.data(), shape, stride, blitz::neverDeleteData);
		TensorPotentialMultiply_Rank1_BandNH(TransformRank, kineticEnergyTensor, -1.0 / (2 * Mass), src, dst);
		
		//blitz::Array<cplx, Rank> kineticEnergyTensor(LaplacianFull.data(), shape, stride, blitz::neverDeleteData);
		//TensorPotentialMultiply_Rank1_Dense(TransformRank, kineticEnergyTensor, -1.0 / (2 * Mass), src, dst);
		
		//MatrixVectorMultiplyBanded(LaplacianBlasBanded, src, dst, -1.0 / (2 * Mass), 0.0);
		
	}
	
	void AdvanceStep(Wavefunction<Rank> &psi, cplx dt)
	{

		blitz::Array<cplx, 1> temp = TempData;
	
		//Scaling factor for matrix-vector multiplication
		cplx scaling = (- I * TimeStep / 2.0) * (- 1.0 / (2 * Mass));

		//Map wavefunction to a rank 3 array, where the middle rank is TransformRank
		blitz::Array<cplx, Rank> fullData = psi.GetData();
		blitz::Array<cplx, 3> data = MapToRank3(fullData, TransformRank, 1);
		
		//Iterate over the array directions which are not propagated by this propagator
		for (int i=0; i<data.extent(0); i++)
		{
			for (int j=0; j<data.extent(2); j++)
			{
				/* 
				 * v is a view of a slice of the wave function
				 * temp is a temporary copy of that slice
				 */
				blitz::Array<cplx, 1> v = data(i, blitz::Range::all(), j);
				temp = v; 
		
				/*
				 * Call BLAS function for fast banded matrix-vector product.
				 *
				 * temp := (1 - i dt / 2 H) psi
				 */
				//TODO:MatrixVectorMultiply for GeneralBanded
				MatrixVectorMultiplyBanded(LaplacianBlasBanded, v, temp, scaling, 1.0);
		
				/*
				 * Then call LAPACK to solve banded system of equations, obtaining 
				 * solution at t+dt/2. We must use a temp array to get the
				 * stride right for LAPACK (stride = 1)
				 *
				 * temp := (1 + i dt / 2 H)^(-1) temp
				 */
				lapack.SolveBandedFactored(BackwardPropagationFactored, BackwardPropagationPivots, temp);
				v = temp;
			}
		}	
	}

	void SetupStep(Wavefunction<Rank> &psi, cplx dt)
	{
		GlobalGrid.reference( psi.GetRepresentation()->GetGlobalGrid(TransformRank).copy() );

		//GlobalGrid = blitz::pow(GlobalGrid+10, 2);

		GlobalGridSize = GlobalGrid.size();
		TimeStep = dt;

		TempData.resize(GlobalGridSize);
		SetupKineticEnergyMatrix();
		SetupPropagationMatrix();
	}

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("difference_order", DifferenceOrder);
		config.Get("transform_rank", TransformRank);
		config.Get("mass", Mass);
	}

	blitz::Array<cplx, 2> GetLaplacianBlasBanded()
	{
		return LaplacianBlasBanded;
	}

	blitz::Array<cplx, 2> GetLaplacianFull()
	{
		return LaplacianFull;
	}

	blitz::Array<cplx, 2> GetBackwardPropagationLapackBanded()
	{
		return BackwardPropagationLapackBanded;
	}

private:
	typedef blitz::linalg::LAPACK<cplx> LAPACK;
	LAPACK lapack;

	int GlobalGridSize;
	blitz::Array<double, 1> GlobalGrid;
	int DifferenceOrder;
	int TransformRank;
	cplx TimeStep;

	double Mass;
	blitz::Array<cplx, 2> LaplacianBlasBanded;  //kinetic energy matrix in the General Banded matrix format
	blitz::Array<cplx, 2> LaplacianFull;  //kinetic energy matrix in full matrix

	blitz::Array<cplx, 2> BackwardPropagationLapackBanded; //Matrix containing (I + i dt H)
	blitz::Array<cplx, 2> BackwardPropagationFactored;     //Matrix containing the LU factorization to the matrix above
	blitz::Array<int, 1> BackwardPropagationPivots;        //Vector containing the pivots to the above factorization

	blitz::Array<cplx, 1> TempData;

	/*
	 * Return the rank 1 representation corresponding to TransformRank
	 */
	CartesianRepresentation<1>::Ptr GetRepresentation(Wavefunction<Rank> &psi)
	{
		typename CombinedRepresentation<Rank>::Ptr combinedRepr = dynamic_pointer_cast< CombinedRepresentation<Rank> >(psi.GetRepresentation());
		return dynamic_pointer_cast< CartesianRepresentation<1> >(combinedRepr->GetRepresentation(TransformRank));
	}


	/* 
	 * Setup second derivative finite difference matrix D^2 of a given order for a equidistant grid.
	 * order must be odd, and > 2.
	 *
	 * reads
	 *	  GlobalGrid, GlobalGridSize, DifferenceOrder
	 * allocates and updates
	 * 	  DiffernceCoefficients, LaplacianHermitianLower
	 *
	 * if order == 5
	 *     / c02  c03  c04   0    0     0   0 \
	 *     | c11  c12  c13  c14   0     0   0 |
	 *     | c20  c21  c22  c23  c24    0   0 |
	 * A = |  0   c30  c31  c32  c33  c34   0 |
	 *     |  0    0   c40  c41  c42  c43  c44 |
	 *     |  0    0    0   c50  c51  c52  c53 |
	 *     \  0    0    0    0   c60  c61  c62 /
	 *
	 * The difference coefficients c are found solving 
	 *           / 0 \  function value
	 *           | 0 |  1st derivative
	 * Bi^T ci = | 1 |  2nd derivative
	 *           | 0 |  ...
	 *           |...|  
	 *           \ 0 /
	 *
	 * Where ci is a row vector i n the above matrix, and
	 *
	 *      /  1  (-3 h_i-3)^1/1! (-3 h_i-3)^2/2! (-3 h_i-3)^3/3! (-3 h_i-3)^4/4! (-3 h_i-3)^5/5! \
	 *      |  1  (-2 h_i-2)^1/1! (-2 h_i-2)^2/2! (-2 h_i-2)^3/3! (-2 h_i-2)^4/4! (-2 h_i-2)^5/5! |
	 *      |  1  (-1 h_i-1)^1/1! (-1 h_i-1)^2/2! (-1 h_i-1)^3/3! (-1 h_i-1)^4/4! (-1 h_i-1)^5/5! |
	 * Bi = |  1         0               0               0               0               0        |
	 *      |  1  ( 1 h_i+1)^1/1! ( 1 h_i+1)^2/2! ( 1 h_i+1)^3/3! ( 1 h_i+1)^4/4! ( 1 h_i+1)^5/5! |
	 *      |  1  ( 2 h_i+2)^1/1! ( 2 h_i+2)^2/2! ( 2 h_i+2)^3/3! ( 2 h_i+2)^4/4! ( 2 h_i+2)^5/5! |
	 *      \  1  ( 3 h_i+3)^1/1! ( 3 h_i+3)^2/2! ( 3 h_i+3)^3/3! ( 3 h_i+3)^4/4! ( 3 h_i+3)^5/5! /
	 *
	 * i.e
	 * B_{i,j} = ((i - (k+1)/2) * h_{i-(k+1)/2} )^j / j!
	 *
	 * where k is the order of the method
	 * and h_{i-l} = x_{i-l} - x_i
	 *
	 * For the simplifying case where the grid is equidistant
	 *
	 */

	/*
	 * Find the difference coefficients c_curIndex, that is, set up a row of the difference matrix
	 */
	blitz::Array<cplx, 1> FindDifferenceCoefficients(int curIndex)
	{
		int k = (DifferenceOrder-1)/2;

		blitz::Array<double, 1> gridDifference(DifferenceOrder);
		gridDifference = 0;

		for (int i=curIndex-k; i<=curIndex+k; i++)
		{
			if (i < 0)
			{
				double endDifference = GlobalGrid(1) - GlobalGrid(0);
				gridDifference(i-curIndex+k) = (endDifference*i + GlobalGrid(0) - GlobalGrid(curIndex));
			}
			else if (i >= GlobalGridSize)
			{
				double endDifference = GlobalGrid(GlobalGridSize-1) - GlobalGrid(GlobalGridSize-2);
				gridDifference(i-curIndex+k) = (endDifference*(i-GlobalGridSize+1) + GlobalGrid(GlobalGridSize-1) - GlobalGrid(curIndex));
			}
			else
			{
				gridDifference(i-curIndex+k) = (GlobalGrid(i) - GlobalGrid(curIndex));
			}
		}

		cout << "curindex = " << curIndex << ", gridDifference = " << ToString(gridDifference) << endl;

		blitz::Array<cplx, 2>  B(DifferenceOrder, DifferenceOrder);
		for (int i=0; i<DifferenceOrder; i++)
		{
			for (int j=0; j<DifferenceOrder; j++)
			{
				B(i, j) = (double)std::pow(gridDifference(i), j) / Factorial(j);
			}
		}

		//Set up input right hand side
		blitz::Array<cplx, 1> differenceCoefficients(DifferenceOrder);
		differenceCoefficients = 0;
		differenceCoefficients(2) = 1;

		blitz::Array<int, 1> pivot(DifferenceOrder);
		lapack.CalculateLUFactorization(B, pivot);
		//LAPACK uses opposite storage model => B is transposed
		lapack.SolveGeneralFactored(LAPACK::TransposeNone, B, pivot, differenceCoefficients);

		//cout << "h = " << ToString(gridDifference) << endl;
		//cout << "c = " << ToString(differenceCoefficients) << endl;
		//cout << endl;

		return differenceCoefficients;
	}


	void SetupKineticEnergyMatrix()
	{
		if (DifferenceOrder <= 2)
		{
			throw std::runtime_error("Can not have 2. derivative finite difference of order < 3");
		}
		if (DifferenceOrder % 2 == 0)
		{
			throw std::runtime_error("Can not have 2. derivative of even order accuracy");
		}

		int k = (DifferenceOrder-1)/2;

		//Set up the difference matrix A
		
		//General Banded BLAS Storage
		LaplacianBlasBanded.resize(GlobalGridSize, DifferenceOrder);
		LaplacianFull.resize(GlobalGridSize, GlobalGridSize);
		LaplacianBlasBanded = 0;
		for (int i=0; i<GlobalGridSize; i++)
		{
			blitz::Array<cplx, 1> differenceCoefficients = FindDifferenceCoefficients(i);

			int startIndex = std::max(0, i-k);
			int endIndex = std::min(GlobalGridSize, i+k+1);
			for (int j=startIndex; j<endIndex; j++)
			{
				int J = j;
				int I = k + i - j;

				LaplacianBlasBanded(J, I) = differenceCoefficients(k + j - i);
				LaplacianFull(i, j) = differenceCoefficients(k + j - i);
			}
		}
	}

	double Factorial(int x)
	{
		return (x < 2) ? (1) : ((double)x * Factorial(x-1));
	}


	/*
	 * Set up the cayley propagation matrices
	 */
	void SetupPropagationMatrix()
	{
		//Number of upper diagonal bands including the main diagonal
		int k = (DifferenceOrder - 1) / 2;

		BackwardPropagationLapackBanded.resize(GlobalGridSize, 3*k + 2);
		BackwardPropagationPivots.resize(GlobalGridSize);
		BackwardPropagationLapackBanded = 0;
		BackwardPropagationPivots = 0;

		//Setup 1 + i dt/2 H
		for (int i = 0; i < GlobalGridSize; i++)
		{
			int startIndex = std::max(0, i-k);
			int endIndex = std::min(GlobalGridSize, i+k+1);
			for (int j=startIndex; j<endIndex; j++)			
			{
				//Upper lapack indices
				int lapackJ = j;
				int lapackI = 2 * k + i - j;

				int blasJ = j; 
				int blasI = k + i - j;

				cplx ham = -1.0 / (2.0 * Mass) * LaplacianBlasBanded(blasJ, blasI);
				
				// 1 + i dt/2 H
				BackwardPropagationLapackBanded(lapackJ, lapackI) = cplx(0.,1.) * TimeStep / 2.0 * ham;
				if (i == j)
				{
					BackwardPropagationLapackBanded(lapackJ, lapackI) += 1;
				}
			}
		}

		//Calculate LU-factorization
		BackwardPropagationFactored.resize(BackwardPropagationLapackBanded.shape());
	 	BackwardPropagationFactored = BackwardPropagationLapackBanded;
	 	lapack.CalculateLUFactorizationBanded(BackwardPropagationFactored, BackwardPropagationPivots);
	}
};

#endif


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
 * Currently only an equidistant grid such as provided by CartesianRepresentation is supported
 *
 * Finite differences are used to approximate the 1-rank kinetic energy operator 1/(2 m) d2/dx2
 *
 */

void TensorPotentialMultiply_Rank1_Band(int rank, blitz::Array<cplx, 1> potential, double scaling, blitz::Array<cplx, 1> &source, blitz::Array<cplx, 1> &dest);
void TensorPotentialMultiply_Rank1_Band(int rank, blitz::Array<cplx, 2> potential, double scaling, blitz::Array<cplx, 2> &source, blitz::Array<cplx, 2> &dest);
void TensorPotentialMultiply_Rank1_Band(int rank, blitz::Array<cplx, 3> potential, double scaling, blitz::Array<cplx, 3> &source, blitz::Array<cplx, 3> &dest);
void TensorPotentialMultiply_Rank1_Band(int rank, blitz::Array<cplx, 4> potential, double scaling, blitz::Array<cplx, 4> &source, blitz::Array<cplx, 4> &dest);

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
		shape(TransformRank) = LaplacianHermitianLower.size();
		blitz::TinyVector<int, Rank> stride = 1;
		for (int j=0; j<TransformRank; j++)
		{
			stride(j) = LaplacianHermitianLower.size();
		}
		blitz::Array<cplx, Rank> kineticEnergyTensor(LaplacianHermitianLower.data(), shape, stride, blitz::neverDeleteData);

		//Use TensorPotential mechanism to perform multiple matrix-vector multiplications
		blitz::Array<cplx, Rank> src = sourcePsi.GetData();
		blitz::Array<cplx, Rank> dst = destPsi.GetData();
		TensorPotentialMultiply_Rank1_Band(TransformRank, kineticEnergyTensor, -1.0 / (2 * Mass), src, dst);
		
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
				MatrixVectorMultiplyHermitianBanded(LaplacianHermitianLower, v, temp, scaling, 1.0);
		
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
		CartesianRepresentation<1>::Ptr repr = this->GetRepresentation(psi);
		GridSpacing = repr->Range(0).Dx;
		GridSize = repr->Range(0).Count;
		TimeStep = dt;

		TempData.resize(GridSize);
		SetupKineticEnergyMatrix();
		SetupPropagationMatrix();
	}

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("difference_order", DifferenceOrder);
		config.Get("transform_rank", TransformRank);
		config.Get("mass", Mass);
	}

	blitz::Array<cplx, 2> GetLaplacianHermitianLower()
	{
		return LaplacianHermitianLower;
	}

	blitz::Array<cplx, 1> GetDifferenceCoefficients()
	{
		return DifferenceCoefficients;
	}

	blitz::Array<cplx, 2> GetBackwardPropagationLapackBanded()
	{
		return BackwardPropagationLapackBanded;
	}

private:
	typedef blitz::linalg::LAPACK<cplx> LAPACK;
	LAPACK lapack;

	int GridSize;
	double GridSpacing;
	int DifferenceOrder;
	int TransformRank;
	cplx TimeStep;

	double Mass;
	blitz::Array<cplx, 1> DifferenceCoefficients;       //coefficients of relative grid indices for d^2/dx^2
	blitz::Array<cplx, 2> LaplacianHermitianLower;  //kinetic energy matrix in the General Banded matrix format

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
	 *	  GridSize, GridSpacing, DifferenceOrder
	 * allocates and updates
	 * 	  DiffernceCoefficients, LaplacianHermitianLower
	 *
	 * if order == 5
	 *     / c2  c3  c4   0   0   0   0 \
	 *     | c1  c2  c3  c4   0   0   0 |
	 *     | c0  c1  c2  c3  c4   0   0 |
	 * A = |  0  c0  c1  c2  c3  c4   0 |
	 *     |  0   0  c0  c1  c2  c3  c4 |
	 *     |  0   0   0  c0  c1  c2  c3 |
	 *     \  0   0   0   0  c0  c1  c2 /
	 *
	 * The difference coefficients c are found solving 
	 *         / 0 \  function value
	 *         | 0 |  1st derivative
	 * B^T c = | 1 |  2nd derivative
	 *         | 0 |  ...
	 *         |...|  
	 *         \ 0 /
	 *
	 * Where 
	 *     /  1  (-3 h)^1/1! (-3 h)^2/2! (-3 h)^3/3! (-3 h)^4/4! (-3 h)^5/5! \
	 *     |  1  (-2 h)^1/1! (-2 h)^2/2! (-2 h)^3/3! (-2 h)^4/4! (-2 h)^5/5! |
	 *     |  1  (-1 h)^1/1! (-1 h)^2/2! (-1 h)^3/3! (-1 h)^4/4! (-1 h)^5/5! |
	 * B = |  1    0           0           0           0           0         |
	 *     |  1  ( 1 h)^1/1! ( 1 h)^2/2! ( 1 h)^3/3! ( 1 h)^4/4! ( 1 h)^5/5! |
	 *     |  1  ( 2 h)^1/1! ( 2 h)^2/2! ( 2 h)^3/3! ( 2 h)^4/4! ( 2 h)^5/5! |
	 *     \  1  ( 3 h)^1/1! ( 3 h)^2/2! ( 3 h)^3/3! ( 3 h)^4/4! ( 3 h)^5/5! /
	 *
	 * i.e
	 * B_{i,j} = ((i - (k+1)/2) * h)^j / j!
	 *
	 * where k is the order of the method
	 *
	 */
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

		//shift of c-indices compared to relative grid point index
		int indexShift = (DifferenceOrder - 1) / 2;

		//Set up coefficient matrix B
		blitz::Array<cplx, 2>  B(DifferenceOrder, DifferenceOrder);
		for (int i=0; i<DifferenceOrder; i++)
		{
			for (int j=0; j<DifferenceOrder; j++)
			{
				B(i, j) = (double)std::pow((i-indexShift)*GridSpacing, j) / Factorial(j);
			}
		}

		//Set up input right hand side
		DifferenceCoefficients.resize(DifferenceOrder);
		DifferenceCoefficients = 0;
		DifferenceCoefficients(2) = 1;

		blitz::Array<int, 1> pivot(DifferenceOrder);
		lapack.CalculateLUFactorization(B, pivot);
		//LAPACK uses opposite storage model => B is transposed
		lapack.SolveGeneralFactored(LAPACK::TransposeNone, B, pivot, DifferenceCoefficients);
		
		//Set up the difference matrix A
		
		//Hermitian Lower Banded Storage
		LaplacianHermitianLower.resize(GridSize, (DifferenceOrder+1)/2);
		LaplacianHermitianLower = 0;
		for (int i=0; i<=indexShift; i++)
		{
			LaplacianHermitianLower(blitz::Range::all(), i) = DifferenceCoefficients(indexShift + i);
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
		int k = (DifferenceOrder + 1) / 2;

		BackwardPropagationLapackBanded.resize(GridSize, 3*k - 1);
		BackwardPropagationPivots.resize(GridSize);
		BackwardPropagationLapackBanded = 0;
		BackwardPropagationPivots = 0;

		//Setup 1 + i dt/2 H
		for (int i = 0; i < GridSize; i++)
		{
			int jMax = std::min(i + k, GridSize);
			for (int j = i; j < jMax; j++)
			{
				//Upper lapack indices
				int Ju = j;
				int Iu = 2 * k - 2 + i - j;

				//Lower lapack indices
				int Jl = i;
				int Il = 2 * k - 2 + j - i;

				//lower hermitian indices
				int Jh = i;
				int Ih = j - i;

				cplx ham = -1.0 / (2.0 * Mass) * LaplacianHermitianLower(Jh, Ih);
				
				//Compute potential matrix element if present
				/*
				if (HasPotential)
				{
					GetPotentialSlice(potentialSlice, i, PotentialVector);
					ham += BSplineObject->BSplineOverlapIntegral(potentialSlice, i, j);
				}
				*/

				// 1 + i dt/2 H
				BackwardPropagationLapackBanded(Ju, Iu) = cplx(0.,1.) * TimeStep / 2.0 * ham;
				BackwardPropagationLapackBanded(Jl, Il) = cplx(0.,1.) * TimeStep / 2.0 * ham;
				if (i == j)
				{
					BackwardPropagationLapackBanded(Jl, Il) += 1;
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


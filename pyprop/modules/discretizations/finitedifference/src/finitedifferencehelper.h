#ifndef FINITEDIFFERENCEHELPER_H
#define FINITEDIFFERENCEHELPER_H

#include <core/common.h>
#include <core/utility/blitzlapack.h>
#include <core/utility/blitzblas.h>
#include <core/utility/blitztricks.h>
#include <core/utility/matrix_conversion.h>


class FiniteDifferenceHelper
{
public:
	typedef blitz::Array<double, 1> GridType;

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

	FiniteDifferenceHelper() {}
	~FiniteDifferenceHelper() {}

	void Setup(GridType grid, int differenceOrder)
	{
		GlobalGrid.reference(grid);
		GlobalGridSize = GlobalGrid.extent(0);
		DifferenceOrder = differenceOrder;
	}

	int GetDifferenceOrder()
	{
		return DifferenceOrder;
	}


	/*
	 * Find the difference coefficients c_curIndex, that is, set up a row of the difference matrix
	 */
	virtual blitz::Array<cplx, 1> FindDifferenceCoefficients(int curIndex)
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

		//cout << "curindex = " << curIndex << ", gridDifference = " << ToString(gridDifference) << endl;

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


	blitz::Array<cplx, 2> SetupLaplacianBlasBanded()
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
	
		//Single Processor
		//General Banded BLAS Storage
		blitz::Array<cplx, 2> laplacianBlasBanded(GlobalGridSize, DifferenceOrder);
		laplacianBlasBanded = 0;
		for (int i=0; i<GlobalGridSize; i++)
		{
			blitz::Array<cplx, 1> differenceCoefficients = FindDifferenceCoefficients(i);

			int startIndex = std::max(0, i-k);
			int endIndex = std::min(GlobalGridSize, i+k+1);
			for (int j=startIndex; j<endIndex; j++)
			{
				int J = j;
				int I = k + i - j;

				laplacianBlasBanded(J, I) = differenceCoefficients(k + j - i);
			}
		}

		return laplacianBlasBanded;
	}

private:
	typedef blitz::linalg::LAPACK<cplx> LAPACK;
	LAPACK lapack;

	GridType GlobalGrid;
	int GlobalGridSize;
	int DifferenceOrder;

	double Factorial(int x)
	{
		return (x < 2) ? (1) : ((double)x * Factorial(x-1));
	}

};


/*
 * Left boundary condition scaling for FD
 *
 * To better account for wavefunction behavior near the origin, the boundary
 * stencils may be changed (effectively accounting for function behavior
 * outside the discretization box). For an order-n rule, there are (n-1)/2 - i
 * free coefficients for stencil (row) i.
 *
 * The free coefficients are scaled according to the given boundaryScaling
 * array, and applied to the grid correspoing grid point inside the grid, i.e.
 * -dr -> dr, -2dr -> 2dr, etc.
 *
 * if order == 5, the free coefficients are as indiciated by 'x' below, and are
 * read from the 1D boundaryScaling array row wise ([0,1,2,...] -> [(0,0),
 * (0,1), (0,2), (1,0), ...])
 *
 *  x x / c02  c03  c04   0    0     0   0  \
 *    x | c11  c12  c13  c14   0     0   0  |
 *      | c20  c21  c22  c23  c24    0   0  |
 *      |  0   c30  c31  c32  c33  c34   0  |
 *      |  0    0   c40  c41  c42  c43  c44 |
 *      |  0    0    0   c50  c51  c52  c53 |
 *      \  0    0    0    0   c60  c61  c62 /
 *
 * Note: if zero boundary conditions are imposed by removing the zero grid
 * point, the offset parameters should be set to -2, otherwise 0.
 *
 * Examples
 * --------
 *  Antisymmetric condition, order 5: M = [-1 -1 -1], offset = 0
 *
 *  Better order-5 rule for Coulomb problems:
 *
 *      D2f(dr) = C*a1*f(dr) + a3*f(dr) + a4*f(2*dr) + a5*f(3*dr)
 *              = (C*a1 + a3)*f(dr) + a4*f(2*dr) + a5*f(3*dr)
 *
 *    c.f. Smyth et. al. 1998. 
 *
 * TODO: Implement for right boundary as well
*/
class FiniteDifferenceHelperCustomBoundary: public FiniteDifferenceHelper
{
public:
	FiniteDifferenceHelperCustomBoundary() {}
	~FiniteDifferenceHelperCustomBoundary() {}

	
	void Setup(GridType grid, int differenceOrder, blitz::Array<double,1> boundaryScaling, int offset)
	{
		FiniteDifferenceHelper::Setup(grid, differenceOrder);
		BoundaryScaling.resize(boundaryScaling.shape());
		BoundaryScaling = boundaryScaling.copy();
		Offset = offset;

		//Check that we got correct number of scaling points
		int k = differenceOrder;
		int numScalingPoints = (k+1) * (k-1) / 8;
		int numInPoints = BoundaryScaling.size();
		if (numInPoints != numScalingPoints)
		{
			std::cout << "Uh-oh, got " << numInPoints << " scaling points, should have been " 
				<< numScalingPoints << endl;
			throw "Got incorrect number of boundary scaling points!";
		}
	}

	/*
	 * Find the difference coefficients c_curIndex, that is, set up a row of the difference matrix
	 */
	virtual blitz::Array<cplx, 1> FindDifferenceCoefficients(int curIndex)
	{
		blitz::Array<cplx, 1> differenceCoefficients =
			FiniteDifferenceHelper::FindDifferenceCoefficients(curIndex);

		int k = this->GetDifferenceOrder();
		int b = (k - 1)/2;

		double m = 0;
		int fdIdx = 0;
		int bIdx = 0;
		int gridIdx = 0;
		for (int j=0; j<(b-curIndex); j++)
		{
			//Calculate boundary scaling index for this FD index
			bIdx = b * curIndex - curIndex * (curIndex - 1)/2 + j;
			m = BoundaryScaling(bIdx);
			
			//Grid and FD stencil index
			fdIdx = j;
			gridIdx = k - 1 + Offset - 2*curIndex - j;
			if (gridIdx < 0)
				continue;

			//cout << "rowIdx = " << curIndex << ", gridIdx = " << (gridIdx) << ", fdIdx = " 
			//	<< fdIdx << ", bIdx = " << bIdx << ", m = " << m << endl;

			//Apply boundary scaling
			differenceCoefficients(gridIdx) += m * differenceCoefficients(fdIdx);
		}

		return differenceCoefficients;
	}

private:
	blitz::Array<double, 1> BoundaryScaling;
	int Offset;
};


#endif


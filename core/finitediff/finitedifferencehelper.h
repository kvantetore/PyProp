#ifndef FINITEDIFFERENCEHELPER_H
#define FINITEDIFFERENCEHELPER_H

#include "../common.h"
#include "../utility/blitzlapack.h"
#include "../utility/blitzblas.h"
#include "../utility/blitztricks.h"


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

blitz::Array<cplx, 2> ConvertMatrixBlasBandedToFull(blitz::Array<cplx, 2> blasBanded)
{
	int N = blasBanded.extent(0);
	int k = (blasBanded.extent(1) - 1) / 2;

	blitz::Array<cplx, 2> full(N, N);
	full = 0;
	for (int i=0; i<N; i++)
	{
		int startIndex = std::max(0, i-k);
		int endIndex = std::min(N, i+k+1);
		for (int j=startIndex; j<endIndex; j++)
		{
			int J = j;
			int I = k + i - j;

			full(i, j) = blasBanded(J, I);
		}
	}

	return full;
}

blitz::Array<cplx, 2> ConvertMatrixBlasBandedToDistributedBanded(blitz::Array<cplx, 2> blasBanded, int localGridSize, int localGridStart)
{
	/* Convert a matrix in the banded blas format to the banded distributed format
	 *
	 *  /  a a   |       \
	 *  |  a a a |       |
	 *  |    a a | a     |
	 *  |      a | a a   |
	 *  |        | a a a | 
	 *  \        |   a a /
	 *
	 */

	int N = blasBanded.extent(0);
	int bandCount = blasBanded.extent(1);
	int k = (bandCount - 1) / 2;

	//Parallel
	blitz::Array<cplx, 2> distributedBanded(localGridSize, bandCount);
	distributedBanded = 0;
	int globalRowStart = std::max(localGridStart-k, 0);
	int globalRowEnd = std::min(localGridStart+localGridSize+k+1, N);
	for (int globalRow=globalRowStart; globalRow<globalRowEnd; globalRow++)
	{
		int globalStartCol = std::max(localGridStart, globalRow-k);
		int globalEndCol = std::min(localGridStart+localGridSize, globalRow+k+1);
		for (int globalCol=globalStartCol; globalCol<globalEndCol; globalCol++)
		{
			int J = globalCol - localGridStart;
			int I = k + globalRow - globalCol;

			int globalJ = globalCol;
			int globalI = I;

			distributedBanded(J, I) = blasBanded(globalJ, globalI);
		}
	}

	return distributedBanded;
};


blitz::Array<cplx, 2> ConvertMatrixBlasBandedToLapackBanded(blitz::Array<cplx, 2> blasBanded)
{
	int N = blasBanded.extent(0);
	int k = (blasBanded.extent(1) - 1) / 2;

	blitz::Array<cplx, 2> lapackBanded(N, 3*k + 2);
	lapackBanded = 0;

	for (int i = 0; i < N; i++)
	{
		int startIndex = std::max(0, i-k);
		int endIndex = std::min(N, i+k+1);
		for (int j=startIndex; j<endIndex; j++)			
		{
			int lapackJ = j;
			int lapackI = 2 * k + i - j;

			int blasJ = j; 
			int blasI = k + i - j;

			lapackBanded(lapackJ, lapackI) = blasBanded(blasJ, blasI);
		}
	}

	return lapackBanded;
}


blitz::Array<cplx, 1> GetDiagonalViewLapackBanded(blitz::Array<cplx, 2> lapackBanded)
{
	int k = (lapackBanded.extent(1) - 2) / 3;
	return lapackBanded(blitz::Range::all(), 2*k);
}


#endif


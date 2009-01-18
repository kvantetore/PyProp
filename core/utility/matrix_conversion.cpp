#include "matrix_conversion.h"

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



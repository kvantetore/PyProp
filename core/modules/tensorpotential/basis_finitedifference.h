#ifndef BASIS_FINITEDIFFERENCE_H
#define BASIS_FINITEDIFFERENCE_H

#include "../common.h"

template<class TBase, int Rank>
void RepresentPotentialInBasisFiniteDifference( Array<TBase, 2> differenceMatrixBandedBlas, Array<TBase, Rank> source, Array<TBase, Rank> dest, Array<int, 2> indexPair, int rank )
{
	typedef Array<TBase, 3> Array3D;
	typedef Array<TBase, 1> Array1D;
	Array3D source3d = MapToRank3(source, rank, 1);
	Array3D dest3d = MapToRank3(dest, rank, 1);

	int preCount = source3d.extent(0);
	int postCount = source3d.extent(2);
	int pairCount = indexPair.extent(0);

	int k = (differenceMatrixBandedBlas.extent(1) - 1) / 2;

	dest3d = 0;
	for (int preIndex=0; preIndex<preCount; preIndex++)
	{
		for (int pairIndex=0; pairIndex<pairCount; pairIndex++)
		{
			int rowIndex = indexPair(pairIndex, 0);
			int colIndex = indexPair(pairIndex, 1);

			if (rowIndex == -1) continue;
			if (colIndex == -1) continue;

			int blasJ = colIndex;
			int blasI = k + rowIndex - colIndex;
			cplx diffCoeff = 0;
			if (0 <= blasI && blasI < differenceMatrixBandedBlas.extent(1))
			{
				diffCoeff = differenceMatrixBandedBlas(blasJ, blasI);
			}

			for (int postIndex=0; postIndex<postCount; postIndex++)
			{
				dest3d(preIndex, pairIndex, postIndex) = diffCoeff * source3d(preIndex, colIndex, postIndex);
			}
		}
	}
}



#endif


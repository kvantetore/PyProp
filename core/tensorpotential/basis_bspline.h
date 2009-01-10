#ifndef BASIS_BSPLINE_H
#define BASIS_BSPLINE_H

#include "../common.h"
#include "../transform/bspline/bspline.h"
#include "../utility/blitztricks.h"

using namespace blitz;

template<class TBase, int Rank>
void RepresentPotentialInBasisBSpline( BSpline::BSpline::Ptr bsplineObject, Array<TBase, Rank> source, Array<TBase, Rank> dest, Array<int, 2> indexPair, std::string storageId, int rank, int differentiation )
{
	typedef Array<TBase, 3> Array3D;
	typedef Array<TBase, 1> Array1D;
	Array3D source3d = MapToRank3(source, rank, 1);
	Array3D dest3d = MapToRank3(dest, rank, 1);

	int preCount = source3d.extent(0);
	int postCount = source3d.extent(2);
	int pairCount = indexPair.extent(0);

	bool isHermitian = storageId == "Herm";
	double scaling = 1;

	for (int preIndex=0; preIndex<preCount; preIndex++)
	{
		for (int pairIndex=0; pairIndex<pairCount; pairIndex++)
		{
			int rowIndex = indexPair(pairIndex, 0);
			int colIndex = indexPair(pairIndex, 1);
			if (isHermitian && rowIndex == colIndex)
			{
				scaling = 0.5;
			}
			else
			{
				scaling = 1.0;
			}

			for (int postIndex=0; postIndex<postCount; postIndex++)
			{
				Array1D sourceSlice = source3d(preIndex, Range::all(), postIndex);
				dest3d(preIndex, pairIndex, postIndex) = scaling * bsplineObject->BSplineGlobalOverlapIntegral(sourceSlice, differentiation, rowIndex, colIndex);
			}
		}
	}
}

#endif


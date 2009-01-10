#ifndef BASIS_REDUCEDSPHERICAL_H
#define BASIS_REDUCEDSPHERICAL_H

#include "../common.h"
#include "../transform/reducedspherical/reducedsphericaltools.h"
#include "../utility/blitztricks.h"

typedef ReducedSpherical::ReducedSphericalTools::Ptr ReducedSphericalToolsPtr;
using namespace blitz;

template<class TBase, int Rank>
void RepresentPotentialInBasisReducedSphericalHarmonic( ReducedSphericalToolsPtr obj, Array<TBase, Rank> source, Array<TBase, Rank> dest, Array<int, 2> indexPair, std::string storageId, int rank, int differentiation )
{
	typedef Array<TBase, 3> Array3D;
	typedef Array<TBase, 1> Array1D;
	Array3D source3d = MapToRank3(source, rank, 1);
	Array3D dest3d = MapToRank3(dest, rank, 1);

	int preCount = source3d.extent(0);
	int thetaCount = source3d.extent(1);
	int postCount = source3d.extent(2);
	int pairCount = indexPair.extent(0);

	
	Array<double, 2> leftAssocLegendre = obj->GetAssociatedLegendrePolynomial();
	Array<double, 2> rightAssocLegendre;
	if (differentiation == 0)
	{
		rightAssocLegendre.reference(obj->GetAssociatedLegendrePolynomial());
	}
	else if (differentiation == 1)
	{
		rightAssocLegendre.reference(obj->GetAssociatedLegendrePolynomialDerivative());
	}
	else 
	{
		throw std::runtime_error("Reduced spherical harmonics does not support higher derivatives than 1");
	}
	Array<double, 1> weights = obj->GetWeights();
	
	dest3d = 0;
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

			for (int thetaIndex=0; thetaIndex<thetaCount; thetaIndex++)
			{
				double leftLegendre = leftAssocLegendre(thetaIndex, rowIndex);
				double rightLegendre = rightAssocLegendre(thetaIndex, colIndex);
				double weight = weights(thetaIndex);
		
				for (int postIndex=0; postIndex<postCount; postIndex++)
				{
					dest3d(preIndex, pairIndex, postIndex) += leftLegendre * source3d(preIndex, thetaIndex, postIndex) * rightLegendre * weight;
				}
			}
		}
	}
}

#endif



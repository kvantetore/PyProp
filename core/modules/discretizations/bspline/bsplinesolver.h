#ifndef BSPLINESOLVER_H
#define BSPLINESOLVER_H

#include <core/common.h>
#include <core/utility/matrix_conversion.h>
#include <core/wavefunction.h>
#include <core/representation/representation.h>
#include "bspline.h"

namespace BSpline 
{

template <int Rank>
class BSplineSolver
{
public:
	typedef blitz::Array<cplx, 2> MatrixType;
	typedef blitz::Array<int, 1> IntVectorType;
	typedef blitz::Array<cplx, Rank> PotentialType;

	blitz::Array<MatrixType, 1> MatrixData;
	PotentialType PotentialData;
	blitz::Array<IntVectorType, 1> PivotData;

	BSplineSolver() {}
	~BSplineSolver() {}

	/* Add potentials
	 *   Potentials should be 
	 *      blas-banded (nonhermitian) in the bspline rank
	 *      diagonal in all other ranks
	 */
	void AddTensorPotential(PotentialType potential)
	{
		if (PotentialData.size() == 0)
		{
			PotentialData.resize(potential.shape());
		}

		PotentialData += potential;
	}

	void Setup(typename Wavefunction<Rank>::Ptr psi, int activeRank, cplx scalingS, cplx scalingH)
	{
		ActiveRank = activeRank;
		typename Representation<Rank>::Ptr repr = psi->GetRepresentation();

		if (PotentialData.size() == 0)
		{
			throw std::runtime_error("No potentials added");
		}

		for (int i=0; i<Rank; i++)
		{
			if (i != ActiveRank)
			{
				if (psi->GetData().extent(i) != PotentialData.extent(i))
				{
					throw std::runtime_error("Invalid potential shape: Potentials must be diagonal in all other ranks than activeRank");
				}
			}
		}

		if (repr->GetDistributedModel()->IsDistributedRank(ActiveRank))
		{
			throw std::runtime_error("ActiveRank can not be distributed");
		}

		//BSpline-basises are not orthogonal
		if (repr->IsOrthogonalBasis(ActiveRank))
		{
			throw std::runtime_error("ActiveRank is an orthogonal basis!");
		}


		//Set up big matrix with correct number of matrices: 
		//  #elements before activeRank * #elements after activeRank
		int numMatrices = psi->GetData().size() / psi->GetData().extent(activeRank);
		MatrixData.resize(numMatrices);
		PivotData.resize(numMatrices);

		//Get Overlap Matrix
		MatrixType overlapBlasBanded = repr->GetGlobalOverlapMatrix(ActiveRank)->GetOverlapBlasBanded();
		MatrixType overlapLapackBanded = ConvertMatrixBlasBandedToLapackBanded(overlapBlasBanded);
		int N = overlapBlasBanded.extent(0);
		int k = (overlapBlasBanded.extent(1) - 1) / 2;

		if (PotentialData.extent(ActiveRank) != overlapBlasBanded.size())
		{
			throw std::runtime_error("Invalid potential shape: Potentials must be blas-banded (nonhermitian) in activeRank");
		}

		//Setup a matrix for all other grid points 
		blitz::Array<cplx, 3> potential3d = MapToRank3(PotentialData, activeRank, 1);
		int matrixIndex = 0;
		int outerCount = potential3d.extent(0);
		int innerCount = potential3d.extent(2);
		for (int outerIndex=0; outerIndex<outerCount; outerIndex++)
		{
			for (int innerIndex=0; innerIndex<innerCount; innerIndex++)
			{
				//Create a view of the current matrix
				PivotData(matrixIndex).resize(N);
				MatrixData(matrixIndex).resize(overlapLapackBanded.shape());
				MatrixType matrix = MatrixData(matrixIndex);
				IntVectorType pivots = PivotData(matrixIndex);

				//Create a view of the current potential
				blitz::Array<cplx, 1> potentialSlice1d = potential3d(outerIndex, blitz::Range::all(), innerIndex);
				blitz::TinyVector<int, 2> potStride((2*k+1) * innerCount, innerCount);
				blitz::TinyVector<int, 2> potShape(overlapBlasBanded.shape());
				MatrixType potentialSliceBlasBanded = MatrixType(potentialSlice1d.data(), potShape, potStride, blitz::neverDeleteData);
				MatrixType potentialLapackBanded = ConvertMatrixBlasBandedToLapackBanded(potentialSliceBlasBanded);

				//Set matrix = (scalingS * S + scalingH + H) 
				//where I is the identity matrix and H is the laplacian plus potentials
				//Set overlap
				matrix = scalingS*overlapLapackBanded + scalingH * potentialLapackBanded;

				//Create LU factorization
			 	lapack.CalculateLUFactorizationBanded(matrix, pivots);	
				
				matrixIndex++;
			}
		}
	}

	void Solve(typename Wavefunction<Rank>::Ptr psi)
	{
		//Map to rank 3 to simplify n-dimensional arrays
		blitz::Array<cplx, 3> psi3d = MapToRank3(psi->GetData(), ActiveRank, 1);

		//Make sure we have a temporary array
		int rankSize = psi->GetData().extent(ActiveRank);
		if (temp.extent(0) != rankSize)
		{
			temp.resize(rankSize);
		}

		//Solve
		int matrixIndex = 0;
		for (int outerIndex=0; outerIndex<psi3d.extent(0); outerIndex++)
		{
			for (int innerIndex=0; innerIndex<psi3d.extent(2); innerIndex++)
			{
				//Create a view of the current matrix
				MatrixType matrix = MatrixData(matrixIndex);
				IntVectorType pivots = PivotData(matrixIndex);
				
				//Create a copy of the current psi to make it unit strided
				temp = psi3d(outerIndex, blitz::Range::all(), innerIndex);

				//Solve for psi
				lapack.SolveBandedFactored(matrix, pivots, temp);

				//copy solution back into psi
				psi3d(outerIndex, blitz::Range::all(), innerIndex) = temp;

				matrixIndex++;
			}
		}
	}


private:
		typedef blitz::linalg::LAPACK<cplx> LAPACK;
		blitz::Array<cplx, 1> temp;
		LAPACK lapack;
		int ActiveRank;
};
	

} //Namespace

#endif


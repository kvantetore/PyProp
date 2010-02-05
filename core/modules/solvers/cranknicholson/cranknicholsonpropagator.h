#ifndef CRANKNICHOLSONPROPAGATOR_H
#define CRANKNICHOLSONPROPAGATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/combinedrepresentation.h"
#include "../representation/cartesianrepresentation.h"
#include "../utility/blitzlapack.h"
#include "../utility/blitzblas.h"
#include "../utility/blitztricks.h"

#include "finitedifferencehelper.h"

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
void TensorPotentialMultiply_Rank1_Distr(int rank, blitz::Array<cplx, 1> potential, double scaling, blitz::Array<cplx, 1> &source, blitz::Array<cplx, 1> &dest, int globalSize, int bands);
void TensorPotentialMultiply_Rank1_Distr(int rank, blitz::Array<cplx, 2> potential, double scaling, blitz::Array<cplx, 2> &source, blitz::Array<cplx, 2> &dest, int globalSize, int bands);
void TensorPotentialMultiply_Rank1_Distr(int rank, blitz::Array<cplx, 3> potential, double scaling, blitz::Array<cplx, 3> &source, blitz::Array<cplx, 3> &dest, int globalSize, int bands);
void TensorPotentialMultiply_Rank1_Distr(int rank, blitz::Array<cplx, 4> potential, double scaling, blitz::Array<cplx, 4> &source, blitz::Array<cplx, 4> &dest, int globalSize, int bands);


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
	CrankNicholsonPropagator() {}
	~CrankNicholsonPropagator() {}

	/* 
	 * Multiplies the kinetic energy operator on sourcePsi and 
	 * adds the result to destPsi.
	 *
	 * Note that it does not multiply the potential
	 */
	void MultiplyKineticEnergyOperator(Wavefunction<Rank> &sourcePsi, Wavefunction<Rank> &destPsi)
	{
		//Check that we're distributed the same way as when SetupStep was called.
		blitz::Array<int, 1> curDistribSrc = sourcePsi.GetRepresentation()->GetDistributedModel()->GetDistribution();
		blitz::Array<int, 1> curDistribDst = destPsi.GetRepresentation()->GetDistributedModel()->GetDistribution();
		bool isDistributed = sourcePsi.GetRepresentation()->GetDistributedModel()->IsDistributedRank(TransformRank);
		if (!blitz::all(Distribution == curDistribSrc && curDistribSrc == curDistribDst))
		{
			cout << "Got distribution " << ToString(curDistribSrc)
			     << ", Expected " << ToString(Distribution) << endl;
			throw std::runtime_error("Invalid distribution encountered in MultiplyKineticEnergyOperator.");
		}

		if (!isDistributed)
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
		else
		{
			//Reshape the kinetic energy matrix into a N-d array suitable for TensorPotentialMultiply
			blitz::TinyVector<int, Rank> shape = 1;
			shape(TransformRank) = LaplacianBlasBanded.size();
			blitz::TinyVector<int, Rank> stride = 1;
			for (int j=0; j<TransformRank; j++)
			{
				stride(j) = LaplacianDistributedBanded.size();
			}
		
			//Use TensorPotential mechanism to perform multiple matrix-vector multiplications
			blitz::Array<cplx, Rank> src = sourcePsi.GetData();
			blitz::Array<cplx, Rank> dst = destPsi.GetData();
		
			int k = (DifferenceOrder-1)/2;
			
			blitz::Array<cplx, Rank> kineticEnergyTensor(LaplacianDistributedBanded.data(), shape, stride, blitz::neverDeleteData);
			TensorPotentialMultiply_Rank1_Distr(TransformRank, kineticEnergyTensor, -1.0 / (2 * Mass), src, dst, GlobalGridSize, k);
		}
		
	}
	
	void AdvanceStep(Wavefunction<Rank> &psi, cplx dt)
	{
		if (psi.GetRepresentation()->GetDistributedModel()->IsDistributedRank(TransformRank))
		{
			throw std::runtime_error("Cannot propagate crank nicholson forward with distributed rank");
		}

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
		//Global Grid
		GlobalGrid.reference( psi.GetRepresentation()->GetGlobalGrid(TransformRank).copy() );
		GlobalGridSize = GlobalGrid.size();
		TimeStep = dt;
		
		//Setup variables for Parallelization
		LocalGrid.reference( psi.GetRepresentation()->GetLocalGrid(TransformRank).copy() );
		LocalGridStart = psi.GetRepresentation()->GetDistributedModel()->GetLocalStartIndex(GlobalGridSize, TransformRank);
		LocalGridSize = LocalGrid.size();
		Distribution.reference(psi.GetRepresentation()->GetDistributedModel()->GetDistribution().copy());

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

	blitz::Array<cplx, 2> GetLaplacianDistributedBanded()
	{
		return LaplacianDistributedBanded;
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

	FiniteDifferenceHelper finiteDifference;

	int GlobalGridSize;
	blitz::Array<double, 1> GlobalGrid;
	int DifferenceOrder;
	int TransformRank;
	cplx TimeStep;

	//Parallelization
	blitz::Array<double, 1> LocalGrid;
	int LocalGridSize;
	int LocalGridStart;
	blitz::Array<int, 1> Distribution;


	double Mass;
	blitz::Array<cplx, 2> LaplacianBlasBanded;  //kinetic energy matrix in the General Banded matrix format
	blitz::Array<cplx, 2> LaplacianDistributedBanded;  //kinetic energy matrix in the distributed banded matrix format
	blitz::Array<cplx, 2> LaplacianFull;  //kinetic energy matrix in full matrix

	blitz::Array<cplx, 2> BackwardPropagationLapackBanded; //Matrix containing (I + i dt H / 2)
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
	 * Setup the laplacian matrix in various formats
	 */
	void SetupKineticEnergyMatrix()
	{
		//Setup the finite difference helper
		finiteDifference.Setup(GlobalGrid, DifferenceOrder);

		//Setup Laplacian Matrix
		LaplacianBlasBanded.reference( finiteDifference.SetupLaplacianBlasBanded() );

		//Convert from Blas Banded format to various other formats
		LaplacianFull.reference( ConvertMatrixBlasBandedToFull(LaplacianBlasBanded) );
		LaplacianDistributedBanded.reference( ConvertMatrixBlasBandedToDistributedBanded(LaplacianBlasBanded, LocalGridSize, LocalGridStart) );
	}

	/*
	 * Set up the cayley propagation matrices
	 */
	void SetupPropagationMatrix()
	{
		//Setup propagation matrix in lapack banded format P = (I + i dt/2 H)
		BackwardPropagationLapackBanded.reference( ConvertMatrixBlasBandedToLapackBanded(LaplacianBlasBanded) );
		//Scale laplacian matrix by i dt/2 and kinetic energy coefficient -0.5/mass
		double kineticCoefficient = -1.0 / (2.0 * Mass);
		BackwardPropagationLapackBanded *= kineticCoefficient * cplx(0.,1.) * TimeStep / 2.0;

		//Add one to diagonal
		blitz::Array<cplx, 1> diagonal = GetDiagonalViewLapackBanded(BackwardPropagationLapackBanded);
		diagonal += 1;

		//Calculate LU-factorization
		BackwardPropagationFactored.resize(BackwardPropagationLapackBanded.shape());
		BackwardPropagationPivots.resize(GlobalGridSize);
	 	BackwardPropagationFactored = BackwardPropagationLapackBanded;
	 	lapack.CalculateLUFactorizationBanded(BackwardPropagationFactored, BackwardPropagationPivots);
	}
};

#endif


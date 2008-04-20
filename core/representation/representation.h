#ifndef REPRESENTATION_H
#define REPRESENTATION_H

#include "../common.h"
#include "../mpi/distributedmodel.h"
#include "overlapmatrix.h"

template<int Rank> class Wavefunction;

template<int Rank>
class Representation
{
public:
	typedef shared_ptr< Representation<Rank> > Ptr;
	typedef typename DistributedModel<Rank>::Ptr DistributedModelPtr;
	typedef boost::shared_ptr< Representation<Rank> > RepresentationPtr;
	
private:
	DistributedModelPtr Distrib;	//Knows how the wavefunction is distributed
	int BaseRank;					//Used when this representation is a part of a combined representation
									//All rank parameters will use the global rank. It is the responsibility
									//of each representation class to translate down to the effective rank
									//(subtract BaseRank)
public:
	//Constructors
	Representation () : BaseRank(0) {}
	virtual ~Representation() {}

	Representation(const Representation<Rank> &other)
	{
		this->Distrib = DistributedModelPtr( new DistributedModel<Rank>(*other.Distrib) );
		this->BaseRank = other.BaseRank;
	}

	inline void SetDistributedModel(DistributedModelPtr distrib)
	{
		Distrib = distrib;
	}
	
	inline DistributedModelPtr GetDistributedModel()
	{
		return Distrib;
	}
	
	blitz::TinyVector<int, Rank> GetInitialShape() 
	{
		blitz::TinyVector<int, Rank> fullShape = GetFullShape();
		return GetDistributedModel()->CreateInitialShape(fullShape);
	}

	int GetBaseRank()
	{
		return BaseRank;
	}

	void SetBaseRank(int baseRank)
	{
		BaseRank = baseRank;
	}

	long int GetId()
	{
		return (long int)this;
	}

	virtual blitz::Array<double, 1> GetLocalGrid(int rank)
	{
		return this->GetDistributedModel()->GetLocalArray(GetGlobalGrid(rank), rank);
	}

	virtual blitz::Array<double, 1> GetLocalWeights(int rank)
	{
		return this->GetDistributedModel()->GetLocalArray(GetGlobalWeights(rank), rank);
	}

	/*
	 * Overlap Matrix related stuff
	 */

	/*
	 * Returns the overlap matrix for a given rank. For orthogonal ranks, this will be a 
	 * diagonal matrix, while for non-orthogonal basises, it will be a banded matrix
	 * (Only compact basises are supported for now)
	 */
	virtual OverlapMatrix::Ptr GetGlobalOverlapMatrix(int rank)
	{
			throw std::runtime_error("OverlapMatrix not implemented for this representation");
	}

	/* 
	 * These are methods from OverlapMatrix, only operated on wavefunctions for all ranks
	 *
	 * Non-Orthogonal basises will typically implement these by calling methods on OverlapMatrix
	 * while orthogonal basises will use the weights to implement it more efficiently
	 */

	virtual void MultiplyOverlap(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi, int rank)
	{
		throw std::runtime_error("MultiplyOverlap not implemented for this representation");
	}

	/*
	 * Multiplies the overlap matrix on the wavefunction in-place Psi := S Psi
	 * S is the tensor product of the overlap matrices for each rank
	 */
	virtual void MultiplyOverlap(Wavefunction<Rank> &psi) 
	{
		throw std::runtime_error("MultiplyOverlap not implemented for this representation");
	}

	/* 
	 * Solves the overlap matrix on the wavefunction in-place S Psi = Psi
	 * S is the tensor product of the overlap matrices for each rank
	 */
	virtual void SolveOverlap(Wavefunction<Rank> &psi) 
	{
		throw std::runtime_error("SolveOverlap not implemented for this representation");
	}

	/* 
	 * Multiplies the sqrt(overlap) or conj(sqrt(S)) matrix on the wavefunction in-place
	 * Psi := sqrt(S) Psi = Psi
	 * S is the tensor product of the overlap matrices for each rank
	 * sqrt(S) = sqrt(S_0) x sqrt(S_1) ...
	 * The square root a positive definite matrix is the cholesky factorization.
	 */
	virtual void MultiplySqrtOverlap(bool conjugate, Wavefunction<Rank> &psi) 
	{
		throw std::runtime_error("MultiplySqrtOverlap not implemented for this representation");
	}

	/* 
	 * Solves the sqrt(overlap) or conj(sqrt(S)) matrix on the wavefunction in-place sqrt(S) Psi = Psi
	 * S is the tensor product of the overlap matrices for each rank
	 * sqrt(S) = sqrt(S_0) x sqrt(S_1) ...
	 * The square root a positive definite matrix is the cholesky factorization.
	 */
	virtual void SolveSqrtOverlap(bool conjugate, Wavefunction<Rank> &psi) 
	{
		throw std::runtime_error("SolveSqrtOverlap not implemented for this representation");
	}

	/*
	 * Checks if the given rank is orthogonal (has diagonal overlap matrix) or not
	 */
	bool IsOrthogonalBasis(int rank)
	{
		return this->GetGlobalOverlapMatrix(rank)->GetSuperDiagonals() == 0;
	}

	//Must override
	virtual blitz::TinyVector<int, Rank> GetFullShape() = 0;
	virtual cplx InnerProduct(const Wavefunction<Rank> &w1, const Wavefunction<Rank> &w2) = 0;
	virtual blitz::Array<double, 1> GetGlobalWeights(int rank) = 0;
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank) = 0;
	virtual void ApplyConfigSection(const ConfigSection &config) = 0;
	virtual RepresentationPtr Copy() = 0;
};

#endif

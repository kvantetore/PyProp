#ifndef REPRESENTATION_H
#define REPRESENTATION_H

#include "../common.h"
#include "../mpi/distributedmodel.h"

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

	/*
	 * Returns the overlap matrix S_{i,j}, where the bands are stored colwise, i.e.
	 * S_{i,j} = overlapMatrix(i, i-j+bw)
	 * where |i - j| <= bw
	 * This means that elements of the same row (of different bands) are stored contiguously in memory
	 */
	virtual blitz::Array<double, 2> GetGlobalOverlapMatrixFullCol(int rank)
	{
		if (IsOrthogonalBasis(rank))
		{
			//For orthogonal basises the overlap matrix is diagonal with the weights on the diagonal
			blitz::Array<double, 1> weights = this->GetLocalWeights(rank);

			blitz::TinyVector<int, 2> shape(weights.extent(0), 1);
			blitz::TinyVector<int, 2> stride(1, 1);
			blitz::Array<double, 2> overlapMatrix(weights.data(), shape, stride, blitz::neverDeleteData);
			return overlapMatrix;
		}
		else
		{
			throw std::runtime_error("OverlapMatrix not implemented for this representation");
		}
	}

	/*
	 * Returns the overlap matrix S_{i,j}, where the bands are stored rowwise, i.e.
	 * S_{i,j} = overlapMatrix(i-j+bw, j)
	 *
	 * This means that elements of the same band is stored contigously in memory
	 */
	virtual blitz::Array<double, 2> GetGlobalOverlapMatrixFullRow(int rank)
	{
		if (IsOrthogonalBasis(rank))
		{
			//For orthogonal basises the overlap matrix is diagonal with the weights on the diagonal
			blitz::Array<double, 1> weights = this->GetLocalWeights(rank);

			blitz::TinyVector<int, 2> shape(1, weights.extent(0));
			blitz::TinyVector<int, 2> stride(1, 1);
			blitz::Array<double, 2> overlapMatrix(weights.data(), shape, stride, blitz::neverDeleteData);
			return overlapMatrix;
		}
		else
		{
			throw std::runtime_error("OverlapMatrix not implemented for this representation");
		}
	}


	virtual int GetOverlapBandwidth(int rank)
	{
		return 1;
	}

	bool IsOrthogonalBasis(int rank)
	{
		return this->GetOverlapBandwidth(rank) == 1;
	}

	//Must override
	virtual blitz::TinyVector<int, Rank> GetFullShape() = 0;
	virtual cplx InnerProduct(const Wavefunction<Rank> &w1, const Wavefunction<Rank> &w2) = 0;
	virtual blitz::Array<double, 1> GetLocalWeights(int rank) = 0;
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank) = 0;
	virtual void ApplyConfigSection(const ConfigSection &config) = 0;
	virtual RepresentationPtr Copy() = 0;
};

#endif

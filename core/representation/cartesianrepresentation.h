#ifndef CARTESIANREPRESENTATION_H
#define CARTESIANREPRESENTATION_H

#include "../common.h"
#include "../mpi/distributedmodel.h"
#include "representation.h"
#include "cartesianrange.h"

template<int Rank>
class CartesianRepresentation : public Representation<Rank>
{
public:
	typedef boost::shared_ptr< CartesianRepresentation<Rank> > Ptr;
	
	blitz::TinyVector<CartesianRange, Rank> Range;
	
	//Constructors
	CartesianRepresentation() {}
	
	CartesianRepresentation(CartesianRange &r0)
	{
		std::cout << "Representation address " << this << std::endl;
	
		for (int i=0; i<Rank; i++)
		{
			Range(i) = r0;
		}
		
		std::cout << "Created cartesian representation of initial shape " << this->GetInitialShape() << std::endl;
	}
	
	CartesianRepresentation(blitz::TinyVector<CartesianRange, Rank> &range)
	{
		for (int i=0; i<Rank; i++)
		{
			Range(i) = range(i);
		}
	}
	
	//Range functions:
	/**
	Range contains the full grid, independent of distributed representation.
	Use GetLocal*() functions to get the portion of the grid local to the current
	processor
	**/
	const CartesianRange& GetRange(int dimension)
	{
		return Range(dimension);
	}
	
	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetLocalGrid(int rank)
	{
		int effectiveRank = rank - this->GetBaseRank();
		return this->GetDistributedModel().GetLocalArray(Range(effectiveRank).GetGrid(), rank);
	}	

	/** 
	Returns the portion of the weights local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetLocalWeights(int rank)
	{
		int effectiveRank = rank - this->GetBaseRank();
		return this->GetDistributedModel().GetLocalArray(Range(effectiveRank).GetWeights(), rank);
	}	

	//Implementation of the Representation interface.
	virtual blitz::TinyVector<int, Rank> GetFullShape();
	virtual cplx InnerProduct(const Wavefunction<Rank> &w1, const Wavefunction<Rank> &w2);
	virtual void ApplyConfigSection(const ConfigSection &cfg);
};

#endif

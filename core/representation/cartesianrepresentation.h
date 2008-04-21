#ifndef CARTESIANREPRESENTATION_H
#define CARTESIANREPRESENTATION_H

#include "../common.h"
#include "../mpi/distributedmodel.h"
#include "cartesianrange.h"
#include "representation.h"

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

	virtual typename Representation<Rank>::RepresentationPtr Copy()
	{
		return typename Representation<Rank>::RepresentationPtr(new CartesianRepresentation<Rank>(*this));
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
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank)
	{
		int effectiveRank = rank - this->GetBaseRank();
		return Range(effectiveRank).GetGrid();
	}	

	/** 
	Returns the portion of the weights local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetGlobalWeights(int rank)
	{
		int effectiveRank = rank - this->GetBaseRank();
		return Range(effectiveRank).GetWeights();
	}	

	/*
	 * Returns the product of all dx
	 */
	double GetScalarWeight()
	{
		double weight = 1;
		for (int i=0; i<Rank; i++)
		{
			weight *= Range(i).Dx;
		}
		return weight;
	}

	//Implementation of the Representation interface.
	virtual blitz::TinyVector<int, Rank> GetFullShape();
	virtual cplx InnerProduct(const Wavefunction<Rank> &w1, const Wavefunction<Rank> &w2);
	virtual void ApplyConfigSection(const ConfigSection &cfg);

	virtual void MultiplyOverlap(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi, int rank);
	virtual void MultiplyOverlap(Wavefunction<Rank> &psi);
	virtual void SolveOverlap(Wavefunction<Rank> &psi);
	virtual void MultiplySqrtOverlap(bool conjugate, Wavefunction<Rank> &psi);
	virtual void SolveSqrtOverlap(bool conjugate, Wavefunction<Rank> &psi);
	virtual OverlapMatrix::Ptr GetGlobalOverlapMatrix(int rank);
};

#endif

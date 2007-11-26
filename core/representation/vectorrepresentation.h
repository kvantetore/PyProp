#ifndef CARTESIANREPRESENTATION_H
#define CARTESIANREPRESENTATION_H

#include "../common.h"
#include "../mpi/distributedmodel.h"
#include "representation.h"
#include "cartesianrange.h"
#include "../utility/blitzblas.h"

class VectorRepresentation : public Representation<1>
{
public:
	typedef boost::shared_ptr< VectorRepresentation > Ptr;
	
	int VectorSize;
	blitz::Array<double, 1> IndexGrid;
	blitz::Array<double, 1> Weights;
	
	//Constructors
	VectorRepresentation() {}
	
	virtual Representation<1>::RepresentationPtr Copy()
	{
		return Representation<1>::RepresentationPtr(new VectorRepresentation(*this));
	}
	
	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank)
	{
		if(IndexGrid.size() == 0) 
		{
			IndexGrid.resize(VectorSize);
			IndexGrid = blitz::tensor::i;
		}

		return IndexGrid;
	}	

	/** 
	Returns the portion of the weights local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetLocalWeights(int rank)
	{
		if (Weights.size() == 0)
		{
			Weights.resize(VectorSize);
			Weights = 1;
		}

		return Weights;
	}	

	//Implementation of the Representation interface.
	virtual blitz::TinyVector<int, 1> GetFullShape()
	{
		return blitz::TinyVector<int, 1>(VectorSize);
	}

	virtual cplx InnerProduct(const Wavefunction<1> &w1, const Wavefunction<1> &w2)
	{
		return VectorInnerProduct(w1.Data, w2.Data);
	}

	virtual void ApplyConfigSection(const ConfigSection &cfg)
	{
		cfg.Get( "vector_size", VectorSize );
	}
};

#endif

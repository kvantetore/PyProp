#ifndef TRANSFORMEDRADIALREPRESENTATION_H
#define TRANSFORMEDRADIALREPRESENTATION_H

#include "../common.h"
#include "transformedrange.h"
#include "representation.h"

/** Represents the wavefunction in a spherical harmonic (l,m) basis.
  * The spherical harmonic of highest order which is represented is 
  * Ylm with l == Range.MaxL, which leaves Range.Count() == (1 + MaxL)**2
  * functions 
  */
class TransformedRadialRepresentation : public Representation<1>
{
public:
	TransformedRange Range;

	//Constructors
	TransformedRadialRepresentation() {}
	virtual ~TransformedRadialRepresentation() {}

	blitz::Array<double, 1> GetFullGrid()
	{
		return Range.GetGrid();
	}
			

	/* ---------- Implementation of Representation<1> interface ----------- */
		
	//Returns the size of the grid 
	virtual blitz::TinyVector<int, 1> GetFullShape()
	{
		return Range.Count;
	}

	/*
	 * Performs an inner product between two radial wavefunctions
	 * This should probably not be called to often, because a faster version
	 * will be available in SphericalRepresentation
	 */
	virtual std::complex<double> InnerProduct(const Wavefunction<1>& w1, const Wavefunction<1>& w2)
	{
		blitz::Range indexRange = this->GetDistributedModel().GetGlobalIndexRange(w1, 0);
		blitz::Array<double, 1> weights(Range.GetWeights());

		cplx innerProd = sum(weights(indexRange) * conj(w1.Data) * w2.Data);
		return innerProd;
	}

	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetLocalGrid(const Wavefunction<1> &psi, int rank)
	{
		blitz::Range indexRange = this->GetDistributedModel().GetGlobalIndexRange(psi, rank);
		return Range.GetGrid()(indexRange);
	}

	/** Apply config, and set up Range
	  */
	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		int N;
		int type;
		TransformedGrid::Parameter param;
		config.Get("n", N);
		config.Get("transform_type", type);
		param.Type = static_cast<TransformedGrid::TransformType>(type);
		config.Get("transform_scaling", param.Scaling);
		
		Range = TransformedRange(param, N);
	}

};

typedef boost::shared_ptr<TransformedRadialRepresentation> TransformedRadialRepresentationPtr;

#endif

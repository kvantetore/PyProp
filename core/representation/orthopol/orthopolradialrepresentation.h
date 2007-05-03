#ifndef ORTHOPOLRADIALREPRESENTATION_H
#define ORTHOPOLRADIALREPRESENTATION_H

#include "../../common.h"
#include "../representation.h"
#include "orthopolrange.h"

class OrthoPolRadialRepresentation : public Representation<1>
{
public:
	OrthoPolRange Range;

	//Constructors
	OrthoPolRadialRepresentation() {}
	virtual ~OrthoPolRadialRepresentation() {}

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
		blitz::Array<double, 1> weights(this->GetDistributedModel().GetLocalArray(Range.GetWeights(), GetBaseRank()));

		cplx innerProd = sum(weights * conj(w1.Data) * w2.Data);
		return innerProd;
	}

	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong transformed radial rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return Range.GetGrid();
	}

	virtual blitz::Array<double, 1> GetLocalWeights(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong transformed radial rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return this->GetDistributedModel().GetLocalArray(Range.GetWeights(), rank);
	}

	/** Apply config, and set up Range
	  */
	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		int N;
		int type;
		OrthoPol::Parameter param;
		config.Get("n", N);
		config.Get("polynomial_type", type);
		param.Type = static_cast<OrthoPol::TransformType>(type);
		config.Get("cutoff", param.Cutoff);
		config.Get("hyperspherical_rank", param.HypersphericalRank);
		
		Range = OrthoPolRange(param, N);
	}

};

typedef boost::shared_ptr<OrthoPolRadialRepresentation> OrthoPolRadialRepresentationPtr;

#endif

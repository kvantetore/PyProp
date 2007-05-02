#ifndef LMREPRESENTATION_H
#define LMREPRESENTATION_H

#include "../../common.h"
#include "../representation.h"
#include "lrange.h"

namespace ReducedSpherical
{

/** Represents the wavefunction in a spherical harmonic (l,m) basis.
  * The spherical harmonic of highest order which is represented is 
  * Ylm with l == Range.MaxL, which leaves Range.Count() == (1 + MaxL)**2
  * functions 
  */
class ReducedSphericalHarmonicRepresentation : public Representation<1>
{
public:
	LRange Range;

	//Constructors
	ReducedSphericalHarmonicRepresentation() {}
	virtual ~ReducedSphericalHarmonicRepresentation() {}

	void SetupRepresentation(int maxL)
	{
		Range = LRange(maxL);
	}

	/* ---------- Implementation of Representation<1> interface ----------- */
		
	//Returns the size of the grid 
	virtual blitz::TinyVector<int, 1> GetFullShape()
	{
		return Range.Count();
	}

	/*
	 * Performs an inner product between two radial wavefunctions
	 * This should probably not be called to often, because a faster version
	 * will be available in SphericalRepresentation
	 */
	virtual std::complex<double> InnerProduct(const Wavefunction<1>& w1, const Wavefunction<1>& w2)
	{
		throw std::runtime_error("Please do not call InnerProduct on the sub representations. Combined representation will take care of proper InnerProduct");
	}

	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetLocalWeights(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong sphharm rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return this->GetDistributedModel().GetLocalArray(Range.GetWeights(), rank);
	}
	
	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong sphharm rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return Range.GetGrid();
	}


	/** Apply config, and set up Range
	  */
	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		int maxl;
		config.Get("maxl", maxl);
		Range = LRange(maxl);
	}

};

typedef boost::shared_ptr<ReducedSphericalHarmonicRepresentation> ReducedSphericalHarmonicRepresentationPtr;

} //Namespace

#endif


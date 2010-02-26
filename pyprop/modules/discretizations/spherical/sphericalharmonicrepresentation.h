#ifndef LMREPRESENTATION_H
#define LMREPRESENTATION_H

#include "../../common.h"
#include "../orthogonalrepresentation.h"
#include "../compressedrepresentation.h"
#include "lmrange.h"

/** Represents the wavefunction in a spherical harmonic (l,m) basis.
  * The spherical harmonic of highest order which is represented is 
  * Ylm with l == Range.MaxL, which leaves Range.Count() == (1 + MaxL)**2
  * functions 
  */
class SphericalHarmonicRepresentation : public CompressedRepresentation
{
public:
	typedef shared_ptr<SphericalHarmonicRepresentation> Ptr;

	LmRange Range;

	//Constructors
	SphericalHarmonicRepresentation() {}
	virtual ~SphericalHarmonicRepresentation() {}

	virtual Representation<1>::RepresentationPtr Copy()
	{
		return Representation<1>::RepresentationPtr(new SphericalHarmonicRepresentation(*this));
	}

	void SetupRepresentation(int maxL)
	{
		Range = LmRange(maxL);
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
		std::cout << "Calculating inner product of spherical harmonics in "
		          << "SphericalHarmonicRepresentation. This should probably be done "
			  << "in the combined representation insted for optimal efficiency"
			  << std::endl;
		
	
		//The inner product in the spherical harmonics representation is
		//simply the sum of the product between each element. This can be shown 
		//from the definition of the expansion coefficients (integration is already
		//taken care of by the spherical harmonics)
		return sum(conj(w1.Data) * w2.Data); 
	}

	/*
	 * Returns Integration weights for the global grid
	 */
	virtual blitz::Array<double, 1> GetGlobalWeights(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong sphharm rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return Range.GetWeights();
	}

	/*
	 * Returns grid points for the global grid. Note that these grid
	 * points are just to satisfy the Representation interface. and
	 * return only the double value of the index number. The actual
	 * l,m values are returned by Get....ExpandedGrid()
	 */
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong sphharm rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return Range.GetIndexGrid();
	}

	/*
	 * Returns the l,m indices for the global grid. See 
	 * CompressedRepresentation for more information
	 */
	virtual blitz::Array<double, 2> GetGlobalExpandedGrid()
	{
		return Range.GetLmGrid();
	}
	
	/** Apply config, and set up Range
	  */
	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		int maxl;
		config.Get("maxl", maxl);
		Range = LmRange(maxl);
	}

};

typedef boost::shared_ptr<SphericalHarmonicRepresentation> SphericalHarmonicRepresentationPtr;

#endif

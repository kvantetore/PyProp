#ifndef SPHERICALHARMONICBASISREPRESENTATION_H
#define SPHERICALHARMONICBASISREPRESENTATION_H

#include "../../common.h"
#include "../orthogonalrepresentation.h"
#include "lmbasisrange.h"

namespace SphericalBasis
{

/** Represents the wavefunction in a spherical harmonic basis.
 *
 * y{l, m}
  */

class SphericalHarmonicBasisRepresentation : public OrthogonalRepresentation
{
public:
	typedef shared_ptr<SphericalHarmonicBasisRepresentation> Ptr;
	
	LmBasisRange Range;

	//Constructors
	SphericalHarmonicBasisRepresentation() {}
	virtual ~SphericalHarmonicBasisRepresentation() {}

	virtual Representation<1>::RepresentationPtr Copy()
	{
		return Representation<1>::RepresentationPtr(new SphericalHarmonicBasisRepresentation(*this));
	}

	/* ---------- Implementation of Representation<1> interface ----------- */
		
	//Returns the size of the grid 
	virtual blitz::TinyVector<int, 1> GetFullShape()
	{
		return Range.Count();
	}

	virtual blitz::Array<double, 1> GetGlobalWeights(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong sphharm rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return Range.GetWeights();
	}
	
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
	virtual void ApplyConfigSection(const ConfigSection &config);

};

typedef boost::shared_ptr<SphericalHarmonicBasisRepresentation> SphericalHarmonicBasisRepresentationPtr;

} //Namespace

#endif


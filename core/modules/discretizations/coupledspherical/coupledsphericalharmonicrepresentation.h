#ifndef COUPLEDSPHERICALHARMONICREPRESENTATION_H
#define COUPLEDSPHERICALHARMONICREPRESENTATION_H

#include "../../common.h"
#include "../orthogonalrepresentation.h"
#include "coupledrange.h"

namespace CoupledSpherical
{

/** Represents the wavefunction in a coupled spherical harmonic basis.
 *
 * y{L,M,l1,l2}
  */

class CoupledSphericalHarmonicRepresentation : public OrthogonalRepresentation
{
public:
	typedef shared_ptr<CoupledSphericalHarmonicRepresentation> Ptr;
	
	CoupledRange Range;

	//Constructors
	CoupledSphericalHarmonicRepresentation() {}
	virtual ~CoupledSphericalHarmonicRepresentation() {}

	virtual Representation<1>::RepresentationPtr Copy()
	{
		return Representation<1>::RepresentationPtr(new CoupledSphericalHarmonicRepresentation(*this));
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

typedef boost::shared_ptr<CoupledSphericalHarmonicRepresentation> CoupledSphericalHarmonicRepresentationPtr;

} //Namespace

#endif


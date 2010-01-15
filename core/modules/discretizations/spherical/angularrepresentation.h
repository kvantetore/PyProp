#ifndef ANGULARREPRESENTATION_H
#define ANGULARREPRESENTATION_H

#include <iostream>
#include "../../common.h"
#include "../orthogonalrepresentation.h"
#include "../compressedrepresentation.h"
#include "omegarange.h"


/** Represents the wavefunction in an angular (theta, phi) basis
  * The distribution of (theta, phi) can be chosen when creating the
  * omega range.
  */

class AngularRepresentation : public CompressedRepresentation
{
public:
	typedef shared_ptr<AngularRepresentation> Ptr;

	OmegaRange Range;                          //The points chosen for the angular grid

	//Constructors:
	AngularRepresentation() {}
	virtual ~AngularRepresentation() {}

	virtual Representation<1>::RepresentationPtr Copy()
	{
		return Representation<1>::RepresentationPtr(new AngularRepresentation(*this));
	}

	void SetupRepresentation(int maxL);

	//Returns the size of the grid
	virtual blitz::TinyVector<int, 1> GetFullShape()
	{
		return blitz::TinyVector<int, 1>(sqr(2*Range.MaxL + 1));
	}

	/*
	 * Performs an inner product between two radial wavefunctions
	 * This should probably not be called to often, because a faster version
	 * will be available in SphericalRepresentation
	 */
	virtual std::complex<double> InnerProduct(const Wavefunction<1>& w1, const Wavefunction<1>& w2)
	{
		throw std::runtime_error("AngularRepresentation::InnerProduct is not implemented");
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
			cout << "Warning: Trying to get the wrong angular rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return Range.GetIndexGrid();
	}

	/** 
	Returns the portion of the weights local to the current processor.
	**/
virtual blitz::Array<double, 1> GetGlobalWeights(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong angular rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return Range.GetWeights();
	}

	/*
	 * Returns the theta, phi values for the global grid. See 
	 * CompressedRepresentation for more information
	 */
	virtual blitz::Array<double, 2> GetGlobalExpandedGrid()
	{
		return Range.GetOmegaGrid();
	}
	
	/** Apply config, and set up Range
	  */
	virtual void ApplyConfigSection(const ConfigSection &config);

};

#endif


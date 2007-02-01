#ifndef ANGULARREPRESENTATION_H
#define ANGULARREPRESENTATION_H

#include <iostream>
#include "../common.h"
#include "representation.h"
#include "omegarange.h"


/** Represents the wavefunction in an angular (theta, phi) basis
  * The distribution of (theta, phi) can be chosen when creating the
  * omega range.
  */
typedef Representation<1> Representation1D;

class AngularRepresentation : public Representation1D
{
public:
	OmegaRange Range;                          //The points chosen for the angular grid
	blitz::Array<double, 1> OmegaWeight;       //The weights used when integrating over the grid

	//Constructors:
	AngularRepresentation() {}
	virtual ~AngularRepresentation() {}

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
		return sum(conj(w1.Data) * w2.Data * OmegaWeight); 
	}

	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetLocalGrid(const Wavefunction<1> &psi, int rank)
	{
		blitz::Range indexRange = this->GetDistributedModel().GetGlobalIndexRange(psi, 0);
		return Range.GetIndexGrid()(indexRange);
	}
	
	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 2> GetLocalOmegaGrid(const Wavefunction<1> &psi)
	{
		//angular representations must be the second rank (rank==1) 
		blitz::Range indexRange = this->GetDistributedModel().GetGlobalIndexRange(psi, 0);
		return Range.GetOmegaGrid()(indexRange, blitz::Range::all());
	}
	
	/** Apply config, and set up Range
	  */
	virtual void ApplyConfigSection(const ConfigSection &config);
};

typedef boost::shared_ptr<AngularRepresentation> AngularRepresentationPtr;

#endif


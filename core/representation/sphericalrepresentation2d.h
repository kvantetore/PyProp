#ifndef SPHERICALREPRESENTATION2D_H
#define SPHERICALREPRESENTATION2D_H

#include "representation.h"
#include "combinedrepresentation.h"
#include "cartesianrange.h"
#include "omegarange.h"
#include "lmrange.h"

/** Represents the wavefunction in an angular (theta, phi) basis
  * The distribution of (theta, phi) can be chosen when creating the
  * omega range.
  */
class AngularRepresentation : public Representation<1>
{
public:
	OmegaRange Range;                          //The points chosen for the angular grid
	blitz::Array<double, 1> OmegaWeight;       //The weights used when integrating over the grid

	//Constructors:
	AngularRepresentation() {}
	virtual ~AngularRepresentation() {}

	//Returns the size of the grid
	virtual blitz::TinyVector<int, 1> GetFullShape()
	{
		return blitz::TinyVector<int, 1>(Range.Count);
	}

	/*
	 * Performs an inner product between two radial wavefunctions
	 * This should probably not be called to often, because a faster version
	 * will be available in SphericalRepresentation
	 */
	cplx InnerProduct(const Wavefunction<1>& w1, const Wavefunction<1>& w2)
	{
		return sum(conj(w1.Data) * w2.Data * OmegaWeight); 
	}

	/** 
	Returns the portion of the grid local to the current processor.
	**/
	blitz::Array<double, 1> GetLocalGrid(const Wavefunction<1> &psi, int rank)
	{
		blitz::Range indexRange = this->GetDistributedModel().GetGlobalIndexRange(psi, rank);
		return Range.GetIndexGrid()(indexRange);
	}	
	
	/** Apply config, and set up Range
	  */
	void ApplyConfigSection(const ConfigSection &config)
	{
		int maxl;
		int gridtype;
		config.Get("maxl", maxl);
		config.Get("angular_type", gridtype);

		//need to determine how many points given maxl and angular type.
		//need to implement the load omega grid functions first
		throw std::runtime_error("AngularRepresentation::ApplyConfigSection is not implemented yet");
	}	
};

/** Represents the wavefunction in a spherical harmonic (l,m) basis.
  * The spherical harmonic of highest order which is represented is 
  * Ylm with l == Range.MaxL, which leaves Range.Count() == (1 + MaxL)**2
  * functions 
  */
class SphericalHarmonicRepresentation : public Representation<1>
{
public:
	LmRange Range;

	//Constructors
	SphericalHarmonicRepresentation() {}
	virtual ~SphericalHarmonicRepresentation() {}
	

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
	cplx InnerProduct(const Wavefunction<1>& w1, const Wavefunction<1>& w2)
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

	/** 
	Returns the portion of the grid local to the current processor.
	**/
	blitz::Array<double, 1> GetLocalGrid(const Wavefunction<1> &psi, int rank)
	{
		blitz::Range indexRange = this->GetDistributedModel().GetGlobalIndexRange(psi, rank);
		return Range.GetIndexGrid()(indexRange);
	}

	/** Apply config, and set up Range
	  */
	void ApplyConfigSection(const ConfigSection &config)
	{
		int maxl;
		config.Get("maxl", maxl);
		Range = LmRange(maxl);
	}

};

/** Specialized CombinedRepresentation for spherical coordinates, so that 
  * we can apply some semantics to the different ranks.
  * Rank == 0 is the radial direction
  * Rank == 1 is the spherical direction.
  * 
  * By implementing the angular and radial directions as completely
  * independent representations, we can more easily change the method of 
  * evauluation in one without affecting the other.
  */
class SphericalRepresentation2D : public CombinedRepresentation<2>
{
public:
	//Constructors
	SphericalRepresentation2D() {}
	virtual ~SphericalRepresentation2D() {}
	

	Representation<1>& GetRadialRepresentation()
	{
		return GetRepresentation(0);
	}

	Representation<1>& GetAngularRepresentation()
	{
		return GetRepresentation(1);
	}

	virtual cplx InnerProduct(const Wavefunction<2>& w1, const Wavefunction<2>& w2)
	{
		throw std::runtime_error("not implemented");
	}

};

#endif


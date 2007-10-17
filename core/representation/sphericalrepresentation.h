#ifndef SPHERICALREPRESENTATION3D_H
#define SPHERICALREPRESENTATION3D_H

#include "../common.h"
#include "representation.h"
#include "combinedrepresentation.h"
#include "cartesianrange.h"


/** Specialized CombinedRepresentation for spherical coordinates, so that 
  * we can apply some semantics to the different ranks.
  *
  * This class has one spherical ranks and some other non-spherical ranks
  * The last rank is the spherical rank, the other ranks can be whatever
  * 
  * By implementing the angular and radial directions as completely
  * independent representations, we can more easily change the method of 
  * evauluation in one without affecting the other.
  */

template<int Rank>
class SphericalRepresentation : public CombinedRepresentation<Rank>
{
public:
	typedef boost::shared_ptr<SphericalRepresentation> Ptr;

	//Constructors
	SphericalRepresentation() {}
	virtual ~SphericalRepresentation() {}

	virtual typename Representation<Rank>::RepresentationPtr Copy()
	{
		return typename Representation<Rank>::RepresentationPtr(new SphericalRepresentation<Rank>(*this));
	}

	Representation1DPtr GetAngularRepresentation()
	{
		return this->GetRepresentation(Rank-1);
	}

	blitz::Array<double, 2> GetLocalAngularGrid();
	virtual std::complex<double> InnerProduct(const Wavefunction<Rank>& w1, const Wavefunction<Rank>& w2);

private:
	std::complex<double> InnerProductSphericalHarmonic(const Wavefunction<Rank>& w1, const Wavefunction<Rank>& w2);
};

#endif


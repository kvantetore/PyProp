#ifndef SPHERICALREPRESENTATION3D_H
#define SPHERICALREPRESENTATION3D_H

#include "../common.h"
#include "representation.h"
#include "combinedrepresentation.h"
#include "cartesianrange.h"


/** Specialized CombinedRepresentation for spherical coordinates, so that 
  * we can apply some semantics to the different ranks.
  * Rank == 0 is the radial direction
  * Rank == 1 is the spherical direction.
  * 
  * By implementing the angular and radial directions as completely
  * independent representations, we can more easily change the method of 
  * evauluation in one without affecting the other.
  */

class SphericalRepresentation3D : public CombinedRepresentation<2>
{
public:
	//Constructors
	SphericalRepresentation3D() {}
	virtual ~SphericalRepresentation3D() {}

	Representation1DPtr GetRadialRepresentation()
	{
		return GetRepresentation(0);
	}

	Representation1DPtr GetAngularRepresentation()
	{
		return GetRepresentation(1);
	}

	blitz::Array<double, 2> GetLocalAngularGrid(const Wavefunction<2>& psi);
	virtual std::complex<double> InnerProduct(const Wavefunction<2>& w1, const Wavefunction<2>& w2);

private:
	std::complex<double> InnerProductSphericalHarmonic(const Wavefunction<2>& w1, const Wavefunction<2>& w2);
	std::complex<double> InnerProductAngular(const Wavefunction<2>& w1, const Wavefunction<2>& w2);
};

#endif


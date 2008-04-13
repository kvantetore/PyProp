#ifndef REDUCEDSPHERICALTRANSFORM_H
#define REDUCEDSPHERICALTRANSFORM_H

#include "../../common.h"
#include "../../wavefunction.h"
#include "../../representation/combinedrepresentation.h"
#include "../../representation/reducedspherical/thetarepresentation.h"
#include "../../representation/reducedspherical/reducedsphericalharmonicrepresentation.h"
#include "reducedsphericaltools.h"

namespace ReducedSpherical
{

/** 
 * Class that performs the transformation between the 
 * angular space and the spherical harmonic space
*/
template<int Rank>
class ReducedSphericalTransform
{
private:
	int BaseRank;

public:
	ReducedSphericalTools transform;

	//Constructors
	ReducedSphericalTransform() :
			BaseRank(-1)
	{}

	void SetupStep(const Wavefunction<Rank> &psi, int baseRank);
	
	/* 
	 * Transforms the wavefunction from Grid Representation to Spherical Harmonic
	 * Representation. 
	 * The first time this function is called, A new data buffer is allocated on the 
	 * wavefunction to accomodate for the out of place transform. 
	 */
	void ForwardTransform(Wavefunction<Rank> &psi);
	void InverseTransform(Wavefunction<Rank> &psi);

	//Create Representations
	ReducedSphericalHarmonicRepresentationPtr CreateSphericalHarmonicRepr();
	ThetaRepresentationPtr CreateAngularRepresentation();

	int GetBaseRank()
	{
		return BaseRank;
	}

	void SetBaseRank(int baseRank)
	{
		BaseRank = baseRank;
	}

};

} //namespace

#endif


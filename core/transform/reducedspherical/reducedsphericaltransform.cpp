#include "reducedsphericaltransform.h"

namespace ReducedSpherical
{

template<int Rank>
void ReducedSphericalTransform<Rank>::SetupStep(const Wavefunction<Rank> &psi, int baseRank)
{
	int sphRank = baseRank;
	SetBaseRank(baseRank);

	// uses the representation of the wavefunction to get MaxL
	typedef CombinedRepresentation<Rank> CmbRepr;
	typedef ReducedSphericalHarmonicRepresentation SphHarmRepr;
	typename CmbRepr::Ptr reprComb = dynamic_pointer_cast< CmbRepr >(psi.GetRepresentation());
	SphHarmRepr::Ptr reprAngular = dynamic_pointer_cast< SphHarmRepr >(reprComb->GetRepresentation(sphRank));
	if (reprAngular == 0) 
	{
		std::cout << "Invalid wavefunction representation, must be ReducedSphericalHarmonicRepresentation" << std::endl;
		throw std::runtime_error("Invalid wavefunction representation");
	}
	transform.Initialize(reprAngular->Range.MaxL, reprAngular->Range.M);
}

/* 
 * Transforms the wavefunction from Grid Representation to Spherical Harmonic
 * Representation. 
 * The first time this function is called, A new data buffer is allocated on the 
 * wavefunction to accomodate for the out of place transform. 
 */
template<int Rank>
void ReducedSphericalTransform<Rank>::ForwardTransform(Wavefunction<Rank> &psi)
{
	blitz::TinyVector<int, Rank> shape = psi.GetData().shape();

	//Get destination buffer
	int sphDataName = psi.GetAvailableDataBufferName(shape);
	if (sphDataName == -1)
	{
		//Allocate dest array if if does not exist
		sphDataName = psi.AllocateData(shape);
	}

	blitz::Array<cplx, Rank> srcData(psi.GetData());
	blitz::Array<cplx, Rank> dstData(psi.GetData(sphDataName));

	//The last rank is the spherical rank
	transform.ForwardTransform(srcData, dstData, GetBaseRank());
	psi.SetActiveBuffer(sphDataName);
}

template<int Rank>
void ReducedSphericalTransform<Rank>::InverseTransform(Wavefunction<Rank> &psi)
{
	blitz::TinyVector<int, Rank> shape = psi.GetData().shape();

	//Get destination buffer
	int angularDataName = psi.GetAvailableDataBufferName(shape);
	if (angularDataName == -1)
	{
		//Allocate dest array if if does not exist
		angularDataName = psi.AllocateData(shape);
	}

	blitz::Array<cplx, Rank> srcData(psi.GetData());
	blitz::Array<cplx, Rank> dstData(psi.GetData(angularDataName));

	//The last rank is the spherical rank
	transform.InverseTransform(srcData, dstData, GetBaseRank());
	psi.SetActiveBuffer(angularDataName);
}

// Creates a spherical harmonic representation based on an Angular representation
template<int Rank>
ReducedSphericalHarmonicRepresentationPtr ReducedSphericalTransform<Rank>::CreateSphericalHarmonicRepr()
{
	ReducedSphericalHarmonicRepresentationPtr repr(new ReducedSphericalHarmonicRepresentation());
	repr->SetupRepresentation( transform.GetLMax() );
	repr->SetBaseRank( GetBaseRank() );
	return repr;
}

template<int Rank>
ThetaRepresentationPtr ReducedSphericalTransform<Rank>::CreateAngularRepresentation()
{
	ThetaRepresentationPtr repr(new ThetaRepresentation());
	repr->SetupRepresentation( transform.GetLMax() );
	repr->SetBaseRank( GetBaseRank() );
	return repr;
}

template class ReducedSphericalTransform<2>;
template class ReducedSphericalTransform<3>;

}


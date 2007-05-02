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
	CmbRepr* reprComb = dynamic_cast< CmbRepr* >(&psi.GetRepresentation());
	SphHarmRepr* reprAngular = dynamic_cast< SphHarmRepr*>(&(*reprComb->GetRepresentation(sphRank)));
	if (reprAngular == 0) 
	{
		std::cout << "Invalid wavefunction representation, must be ReducedSphericalHarmonicRepresentation" << std::endl;
		throw std::runtime_error("Invalid wavefunction representation");
	}
	transform.Initialize(reprAngular->Range.MaxL);
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
	if (SphDataName == -1)
	{
		//This is the first time this function is called. That means the
		//active data buffer is the grid representation and we must allocate
		//a new buffer for spherical harmonic representation
		blitz::TinyVector<int, Rank> shape = psi.GetData().shape();
		
		AngularDataName = psi.GetActiveBufferName();
		SphDataName = psi.AllocateData(shape);
	}

	//Check that we've got the right active databuffer
	if (AngularDataName != psi.GetActiveBufferName())
	{
		throw std::runtime_error("Active databuffer is not what is should be...");
	}

	blitz::Array<cplx, Rank> srcData(psi.GetData());
	blitz::Array<cplx, Rank> dstData(psi.GetData(SphDataName));

	//The last rank is the spherical rank
	transform.ForwardTransform(srcData, dstData, GetBaseRank());
	psi.SetActiveBuffer(SphDataName);
}

template<int Rank>
void ReducedSphericalTransform<Rank>::InverseTransform(Wavefunction<Rank> &psi)
{
	if (AngularDataName == -1)
	{
		//This is the first time any transform function is called on this
		//class. That means the active data buffer is the spherical harmonic
		//representation, and we must allocate a data buffer for the
		//grid data.
		blitz::TinyVector<int, Rank> shape = psi.GetData().shape();
		
		AngularDataName = psi.AllocateData(shape);
		SphDataName = psi.GetActiveBufferName();
	}

	if (SphDataName != psi.GetActiveBufferName())
	{
		throw std::runtime_error("Active data buffer is not the lm buffer");
	}
	
	blitz::Array<cplx, Rank> srcData(psi.GetData());
	blitz::Array<cplx, Rank> dstData(psi.GetData(AngularDataName));

	//The last rank is the spherical rank
	transform.InverseTransform(srcData, dstData, GetBaseRank());
	psi.SetActiveBuffer(AngularDataName);
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


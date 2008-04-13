#ifndef SPHERICALTRANSFORM_H
#define SPHERICALTRANSFORM_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/angularrepresentation.h"
#include "../representation/combinedrepresentation.h"
#include "../representation/sphericalharmonicrepresentation.h"
#include "shtools.h"


/** class that performs the transformation between the angular space and the spherical harmonic space
*/
template<int Rank>
class SphericalTransform
{
public:
	SphericalTransformTensorGrid transform;

	//Constructors
	SphericalTransform() {}

	// get the legendre polinomials evaluated at cos(theta) in shtools object
	// alternative: use the constructor with the Wavefunction argument
	void SetupStep(const Wavefunction<Rank> &psi)
	{
		throw std::runtime_error("SphericalTransform temporarily out of service...");
		/*
		// uses the representation of the wavefunction to get MaxL
		typedef CombinedRepresentation<Rank> CombRepr;
		typedef SphericalHarmonicRepresentation SphHarmRepr;
		typename CombRepr::Ptr reprSphere = dynamic_pointer_cast<CombRepr>(psi.GetRepresentation());
		SphHarmRepr::Ptr reprAngular = dynamic_pointer_cast<SphHarmRepr>(reprSphere->GetAngularRepresentation());
		if (reprAngular == 0) 
		{
			std::cout << "Invalid wavefunction representation, must be SphericalHarmonicRepresentation" << std::endl;
			throw std::runtime_error("Invalid wavefunction representation");
		}
		transform.Initialize(reprAngular->Range.MaxL);
		*/
	}
	
	/* 
	 * Transforms the wavefunction from Grid Representation to Spherical Harmonic
	 * Representation. 
	 * The first time this function is called, A new data buffer is allocated on the 
	 * wavefunction to accomodate for the out of place transform. 
	 */
	void ForwardTransform(Wavefunction<Rank> &psi)
	{
		//Get shape of dest array
		blitz::TinyVector<int, Rank> shape = psi.GetData().shape();
		int lmSize = transform.GetAssociatedLegendrePolynomial().extent(1);
		shape(Rank-1) = lmSize;

		//Get destination buffer
		int lmDataName = psi.GetAvailableDataBufferName(shape);
		if (lmDataName == -1)
		{
			//Allocate dest array if if does not exist
			lmDataName = psi.AllocateData(shape);
		}

		blitz::Array<cplx, Rank> srcData(psi.GetData());
		blitz::Array<cplx, Rank> dstData(psi.GetData(lmDataName));

		//The last rank is the spherical rank
		transform.ForwardTransform(srcData, dstData, Rank-1);
		psi.SetActiveBuffer(lmDataName);
	}

	void InverseTransform(Wavefunction<Rank> &psi)
	{
		//Get shape of dest array
		blitz::TinyVector<int, Rank> shape = psi.GetData().shape();
		int omegaSize = transform.GetOmegaGrid().extent(0);
		shape(Rank-1) = omegaSize;

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
		transform.InverseTransform(srcData, dstData, Rank-1);
		psi.SetActiveBuffer(angularDataName);
	}

	// Creates a spherical harmonic representation based on an Angular representation
	SphericalHarmonicRepresentation::Ptr CreateSphericalHarmonicRepr()
	{
		SphericalHarmonicRepresentation::Ptr repr(new SphericalHarmonicRepresentation());
		repr->SetupRepresentation( transform.GetLMax() );
		repr->SetBaseRank( Rank - 1 );
		return repr;
	}

	AngularRepresentation::Ptr CreateAngularRepresentation()
	{
		AngularRepresentation::Ptr repr(new AngularRepresentation());
		repr->SetupRepresentation( transform.GetLMax() );
		repr->SetBaseRank( Rank - 1 );
		return repr;
	}

};

#endif


#ifndef SPHERICALTRANSFORM_H
#define SPHERICALTRANSFORM_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/sphericalrepresentation.h"
#include "../representation/angularrepresentation.h"
#include "../representation/sphericalharmonicrepresentation.h"
#include "shtools.h"


/** class that performs the transformation between the angular space and the spherical harmonic space
*/
template<int Rank>
class SphericalTransform
{
private:
	int LmDataName;			//
	int AngularDataName;

public:
	SphericalTransformTensorGrid transform;

	//Constructors
	SphericalTransform() :
			LmDataName(-1),
			AngularDataName(-1)
	{}

	// get the legendre polinomials evaluated at cos(theta) in shtools object
	// alternative: use the constructor with the Wavefunction argument
	void SetupStep(const Wavefunction<Rank> &psi)
	{
		// uses the representation of the wavefunction to get MaxL
		typedef SphericalRepresentation<Rank> SphRepr;
		typedef SphericalHarmonicRepresentation SphHarmRepr;
		typename SphRepr::Ptr reprSphere = dynamic_pointer_cast<SphRepr>(psi.GetRepresentation());
		SphHarmRepr::Ptr reprAngular = dynamic_pointer_cast<SphHarmRepr>(reprSphere->GetAngularRepresentation());
		if (reprAngular == 0) 
		{
			std::cout << "Invalid wavefunction representation, must be SphericalHarmonicRepresentation" << std::endl;
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
	void ForwardTransform(Wavefunction<Rank> &psi)
	{
		if (LmDataName == -1)
		{
			//This is the first time this function is called. That means the
			//active data buffer is the grid representation and we must allocate
			//a new buffer for spherical harmonic representation
			AngularDataName = psi.GetActiveBufferName();

			blitz::TinyVector<int, Rank> shape = psi.GetData().shape();
			int lmSize = transform.GetAssociatedLegendrePolynomial().extent(1);
			shape(Rank-1) = lmSize;
			LmDataName = psi.AllocateData(shape);
		}

		//Check that we've got the right active databuffer
		if (AngularDataName != psi.GetActiveBufferName())
		{
			throw std::runtime_error("Active databuffer is not what is should be...");
		}
	
		blitz::Array<cplx, Rank> srcData(psi.GetData());
		blitz::Array<cplx, Rank> dstData(psi.GetData(LmDataName));

		//The last rank is the spherical rank
		transform.ForwardTransform(srcData, dstData, Rank-1);
		psi.SetActiveBuffer(LmDataName);
	}

	void InverseTransform(Wavefunction<Rank> &psi)
	{
		if (AngularDataName == -1)
		{
			//This is the first time any transform function is called on this
			//class. That means the active data buffer is the spherical harmonic
			//representation, and we must allocate a data buffer for the
			//grid data.
			blitz::TinyVector<int, Rank> shape = psi.GetData().shape();
			int omegaSize = transform.GetOmegaGrid().extent(0);
			shape(Rank-1) = omegaSize;
	
			AngularDataName = psi.AllocateData(shape);
			LmDataName = psi.GetActiveBufferName();
		}

		if (LmDataName != psi.GetActiveBufferName())
		{
			throw std::runtime_error("Active data buffer is not the lm buffer");
		}
		
		blitz::Array<cplx, Rank> srcData(psi.GetData());
		blitz::Array<cplx, Rank> dstData(psi.GetData(AngularDataName));
	
		//The last rank is the spherical rank
		transform.InverseTransform(srcData, dstData, Rank-1);
		psi.SetActiveBuffer(AngularDataName);
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


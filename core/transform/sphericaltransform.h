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
	void SetupStep(const Wavefunction<2> &psi)
	{
		// uses the representation of the wavefunction to get MaxL
		typedef SphericalRepresentation<Rank> SphRepr;
		typedef SphericalHarmonicRepresentation SphHarmRepr;
		SphRepr* reprSphere = dynamic_cast< SphRepr* >(&psi.GetRepresentation());
		SphHarmRepr* reprAngular = dynamic_cast<SphHarmRepr*>(&(*reprSphere->GetAngularRepresentation()));
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
	void ForwardTransform(Wavefunction<2> &psi)
	{
		if (LmDataName == -1)
		{
			//This is the first time this function is called. That means the
			//active data buffer is the grid representation and we must allocate
			//a new buffer for spherical harmonic representation
			AngularDataName = psi.GetActiveBufferName();

			int rSize = psi.GetData().extent(0);
			int lmSize = transform.GetAssociatedLegendrePolynomial().extent(1);
			blitz::TinyVector<int, 2> shape(rSize, lmSize);
			LmDataName = psi.AllocateData(shape);
		}

		//Check that we've got the right active databuffer
		if (AngularDataName != psi.GetActiveBufferName())
		{
			throw std::runtime_error("Active databuffer is not what is should be...");
		}
	
		blitz::Array<cplx, 2> srcData(psi.GetData());
		blitz::Array<cplx, 2> dstData(psi.GetData(LmDataName));

		//The last rank is the spherical rank
		transform.ForwardTransform(srcData, dstData, Rank-1);
		psi.SetActiveBuffer(LmDataName);
	}

	void InverseTransform(Wavefunction<2> &psi)
	{
		if (AngularDataName == -1)
		{
			//This is the first time any transform function is called on this
			//class. That means the active data buffer is the spherical harmonic
			//representation, and we must allocate a data buffer for the
			//grid data.
			int rSize = psi.GetData().extent(0);
			int omegaSize = transform.GetOmegaGrid().extent(0);
			blitz::TinyVector<int, 2> shape(rSize, omegaSize);

			AngularDataName = psi.AllocateData(shape);
			LmDataName = psi.GetActiveBufferName();
		}

		if (LmDataName != psi.GetActiveBufferName())
		{
			throw std::runtime_error("Active data buffer is not the lm buffer");
		}
		
		blitz::Array<cplx, 2> srcData(psi.GetData());
		blitz::Array<cplx, 2> dstData(psi.GetData(AngularDataName));
	
		//The last rank is the spherical rank
		transform.InverseTransform(srcData, dstData, Rank-1);
		psi.SetActiveBuffer(AngularDataName);
	}

	// Creates a spherical harmonic representation based on an Angular representation
	SphericalHarmonicRepresentationPtr CreateSphericalHarmonicRepr()
	{
		SphericalHarmonicRepresentationPtr repr(new SphericalHarmonicRepresentation());
		repr->SetupRepresentation( transform.GetLMax() );
		return repr;
	}

	boost::shared_ptr<AngularRepresentation> CreateAngularRepresentation()
	{
		AngularRepresentationPtr repr(new AngularRepresentation());
		repr->SetupRepresentation( transform.GetLMax() );
		return repr;
	}

};

#endif


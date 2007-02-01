#include "sphericalrepresentation3d.h"

#include "angularrepresentation.h"
#include "sphericalharmonicrepresentation.h"
#include "cartesianrepresentation.h"
#include "transformedradialrepresentation.h"

typedef CartesianRepresentation<1>::Ptr CartesianRepresentation1DPtr;
using boost::dynamic_pointer_cast;
using namespace blitz;

Array<double, 2> SphericalRepresentation3D::GetLocalAngularGrid(const Wavefunction<2>& psi)
	{
	//Create sliced psi
	TinyVector<int, 1> pos = 0;
	Wavefunction<1> *slicedPsi = CreateSlicedWavefunction(psi, 1, pos);
	Array<double, 2> grid;
	
	AngularRepresentationPtr angRepr = dynamic_pointer_cast<AngularRepresentation>(GetAngularRepresentation());
	if (angRepr)
	{
		grid.reference(angRepr->GetLocalOmegaGrid(*slicedPsi));
	}
	else
	{
		SphericalHarmonicRepresentationPtr sphRepr = dynamic_pointer_cast<SphericalHarmonicRepresentation>(GetAngularRepresentation()); 
		if (sphRepr)
		{
			grid.reference(sphRepr->GetLocalLmGrid(*slicedPsi));
		}
	}
	if (grid.size() == 0)
	{
		throw std::runtime_error("Invalid Angular Representation for SphericalDynamicPotentialEvaluator");
	}
	
	delete slicedPsi;
	return grid;
}

cplx SphericalRepresentation3D::InnerProduct(const Wavefunction<2>& w1, const Wavefunction<2>& w2)
{
	if (dynamic_pointer_cast<AngularRepresentation>(GetAngularRepresentation()))
	{
		return InnerProductAngular(w1, w2);
	}
	if (dynamic_pointer_cast<SphericalHarmonicRepresentation>(GetAngularRepresentation()))
	{
		return InnerProductSphericalHarmonic(w1, w2);
	}
	
	throw std::runtime_error("Spherical inner product is only implemented for SphericalHarmonic and Angular representation");
}



cplx SphericalRepresentation3D::InnerProductSphericalHarmonic(const Wavefunction<2>& w1, const Wavefunction<2>& w2)
{
	Array<cplx, 2> d1(w1.GetData());
	Array<cplx, 2> d2(w2.GetData());

	//Get representation
	Representation1DPtr repr = GetRadialRepresentation();
	CartesianRepresentation1DPtr cartRepr = dynamic_pointer_cast< CartesianRepresentation<1> >(repr);
	TransformedRadialRepresentationPtr transRepr = dynamic_pointer_cast< TransformedRadialRepresentation >(repr);
	if (cartRepr)
	{
		double weight = cartRepr->GetRange(0).Dx;

		cplx prod =  sum(conj(d1) * d2);
		return prod * weight;
	} 
	else if (transRepr)
	{
		using tensor::i;
		using tensor::j;

		Array<double, 1> weights = transRepr->Range.GetWeights();
		cplx prod = sum(weights(i) * conj(d1(i,j)) * d2(i,j));
		return prod;
	}
	else
	{
		throw std::runtime_error("Invalid Radial Representation");
	}
}


cplx SphericalRepresentation3D::InnerProductAngular(const Wavefunction<2>& w1, const Wavefunction<2>& w2)
{
	//TODO IMPLEMENT CORRECTLY
	std::cout << "Warning Angular inner product doesn't work properly. "
              << "Consider fixing it or use the inner product of spherical "
			  << "harmonic representation instead" << std::endl;
	using namespace blitz;

	CartesianRepresentation<1> *radialRepr = dynamic_cast< CartesianRepresentation<1>* >(&(*GetRadialRepresentation()));
	AngularRepresentation *angRepr = dynamic_cast<AngularRepresentation*>(&(*GetAngularRepresentation()));

	if (radialRepr == 0 || angRepr == 0)
	{
		throw std::runtime_error("Inner product is only implemented when both representations are in grid-space");
	}

	int thetaCount = 2 * angRepr->Range.MaxL + 1;
	int phiCount = 2 * angRepr->Range.MaxL + 1;
	double dphi = 2 * M_PI / (double)phiCount;
	
	const Array<cplx, 2> d1(w1.GetData());
	const Array<cplx, 2> d2(w2.GetData());
	const Array<double, 1> w(angRepr->Range.GetWeights());
			
	cplx innerprod = 0.0;
	cplx innerPhi = 0.0;
	cplx weight = 0.0;
	cplx prod;
	for (int i=0;i<w1.GetData().extent(0);i++)
	{
		int omegaIndex = 0;
		for(int thetaIndex=0; thetaIndex<thetaCount; thetaIndex++)
		{
			weight = w(thetaIndex);
			innerPhi = 0.0;
			for(int phiIndex=0; phiIndex<phiCount; phiIndex++)
			{	
				prod = conj(d1(i, omegaIndex)) * d2(i, omegaIndex);
				innerPhi += prod;
				omegaIndex++;
			}
			innerprod = innerPhi * weight;
		}
	}

	//scale by dr
	innerprod *= radialRepr->GetRange(0).Dx;
	//scale by dphi
	innerprod *= dphi;

	return innerprod;	
}

	



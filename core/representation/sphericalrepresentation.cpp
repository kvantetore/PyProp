#include "sphericalrepresentation.h"

#include "angularrepresentation.h"
#include "sphericalharmonicrepresentation.h"

using boost::dynamic_pointer_cast;

template<int Rank>
blitz::Array<double, 2> SphericalRepresentation<Rank>::GetLocalAngularGrid()
{
	blitz::Array<double, 2> grid;
	AngularRepresentationPtr angRepr = dynamic_pointer_cast<AngularRepresentation>(GetAngularRepresentation());
	if (angRepr)
	{
		grid.reference(angRepr->GetLocalOmegaGrid());
	}
	else
	{
		SphericalHarmonicRepresentationPtr sphRepr = dynamic_pointer_cast<SphericalHarmonicRepresentation>(GetAngularRepresentation()); 
		if (sphRepr)
		{
			grid.reference(sphRepr->GetLocalLmGrid());
		}
	}
	if (grid.size() == 0)
	{
		throw std::runtime_error("Invalid Angular Representation for SphericalDynamicPotentialEvaluator");
	}
	return grid;
}

template<int Rank>
cplx SphericalRepresentation<Rank>::InnerProduct(const Wavefunction<Rank>& w1, const Wavefunction<Rank>& w2)
{
	if (dynamic_pointer_cast<AngularRepresentation>(GetAngularRepresentation()))
	{
		throw std::runtime_error("Inner product not supported in angular representation, please transform to spherical harmonic representatio");
	}
	if (dynamic_pointer_cast<SphericalHarmonicRepresentation>(GetAngularRepresentation()))
	{
		return InnerProductSphericalHarmonic(w1, w2);
	}
	
	throw std::runtime_error("Spherical inner product is only implemented for SphericalHarmonic representation");
}


template<int Rank>
cplx SphericalRepresentation<Rank>::InnerProductSphericalHarmonic(const Wavefunction<Rank>& w1, const Wavefunction<Rank>& w2)
{
	blitz::Array<cplx, Rank> d1(w1.GetData());
	blitz::Array<cplx, Rank> d2(w2.GetData());

	blitz::TinyVector< blitz::Array<double, 1>, Rank-1> weights;
	for (int i=0; i<Rank-1; i++)
	{
		weights(i).reference(CombinedRepresentation<Rank>::GetLocalWeights(i));
	}

	double weight = 1;
	cplx innerProduct = 0;
	
	typename blitz::Array<cplx, Rank>::iterator it1 = d1.begin();
	typename blitz::Array<cplx, Rank>::iterator it2 = d2.begin();
	for (int linearCount=0; linearCount<w1.Data.size(); linearCount++)
	{
		weight = 1;
		for (int curRank=0; curRank<Rank-1; curRank++)
		{
			weight *= weights(curRank)(it1.position()(curRank));
		}
		
		innerProduct += conj(*it1) * (*it2) * weight;
		it1++;
		it2++;
	}

	return innerProduct;
}

template class SphericalRepresentation<2>;
template class SphericalRepresentation<3>;


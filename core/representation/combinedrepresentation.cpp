#include "combinedrepresentation.h"
#include "cartesianrepresentation.h"

/*----------------------------------------------------------------------------
                Implementation of the Representation interface
  ----------------------------------------------------------------------------*/
template<int Rank> blitz::Array<double, 1> CombinedRepresentation<Rank>::
GetGlobalGrid(int rank)
{
	//Get the local grid in the specified rank
	return GetRepresentation(rank)->GetGlobalGrid(rank);
}

template<int Rank> blitz::Array<double, 1> CombinedRepresentation<Rank>::
GetLocalWeights(int rank)
{
	//Get the local grid in the specified rank
	return GetRepresentation(rank)->GetLocalWeights(rank);
}

template<int Rank> blitz::TinyVector<int, Rank> CombinedRepresentation<Rank>
::GetFullShape() 
{
	blitz::TinyVector<int, Rank> shape;
	for (int i=0;i<Rank;i++)
	{
		shape(i) = GetRepresentation(i)->GetFullShape()(0);
	}
	return shape;
}

//Get the representation of the specified rank
template<int Rank> Representation1DPtr CombinedRepresentation<Rank>
::GetRepresentation(int rank)
{
	return Representations[rank];
}

//Change the representation of the specified rank
template<int Rank> void CombinedRepresentation<Rank>
::SetRepresentation(int rank, Representation1DPtr repr)
{
	Representations[rank] = repr;
}

template<int Rank> void CombinedRepresentation<Rank>
::ApplyConfigSection(const ConfigSection &config) 
{
	//Do this manually for each sub-representation
}

template<int Rank> cplx CombinedRepresentation<Rank>
::InnerProduct(const Wavefunction<Rank>& w1, const Wavefunction<Rank>& w2)
{
	blitz::Array<cplx, Rank> d1(w1.GetData());
	blitz::Array<cplx, Rank> d2(w2.GetData());

	blitz::TinyVector< blitz::Array<double, 1>, Rank> weights;
	for (int i=0; i<Rank; i++)
	{
		weights(i).reference(CombinedRepresentation<Rank>::GetLocalWeights(i));
	}

	double weight = 1;
	cplx innerProduct = 0;
	
	typename blitz::Array<cplx, Rank>::iterator it1 = d1.begin();
	typename blitz::Array<cplx, Rank>::iterator it2 = d2.begin();
	for (int linearCount=0; linearCount<w1.Data.size(); linearCount++)
	{
		weight = 1.;
		for (int curRank=0; curRank<Rank; curRank++)
		{
			weight *= weights(curRank)(it1.position()(curRank));
		}
		
		innerProduct += conj(*it1) * (*it2) * weight;
		it1++;
		it2++;
	}

	return innerProduct;
}



template class CombinedRepresentation<1>;
template class CombinedRepresentation<2>;
template class CombinedRepresentation<3>;


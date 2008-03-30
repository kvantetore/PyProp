#include "combinedrepresentation.h"
#include "cartesianrepresentation.h"

#include "../utility/blitzblas.h"

template<int Rank>
CombinedRepresentation<Rank>::CombinedRepresentation() 
{
	Algorithm = 1;
}

template<int Rank>
CombinedRepresentation<Rank>::~CombinedRepresentation() 
{
}


template<int Rank>
typename Representation<Rank>::RepresentationPtr CombinedRepresentation<Rank>::Copy()
{
	return typename Representation<Rank>::RepresentationPtr(new CombinedRepresentation<Rank>(*this));
}


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

template<int Rank> blitz::Array<double, 2> CombinedRepresentation<Rank>::
GetGlobalOverlapMatrixFullRow(int rank)
{
	return GetRepresentation(rank)->GetGlobalOverlapMatrixFullRow(rank);
}

template<int Rank> blitz::Array<double, 2> CombinedRepresentation<Rank>::
GetGlobalOverlapMatrixFullCol(int rank)
{
	return GetRepresentation(rank)->GetGlobalOverlapMatrixFullCol(rank);
}

template<int Rank> int CombinedRepresentation<Rank>::
GetOverlapBandwidth(int rank)
{
	return GetRepresentation(rank)->GetOverlapBandwidth(rank);
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
	//return Representations[rank];
	return Representations(rank);
}

//Change the representation of the specified rank
template<int Rank> void CombinedRepresentation<Rank>
::SetRepresentation(int rank, Representation1DPtr repr)
{
	Representations(rank) = repr;
}

template<int Rank> void CombinedRepresentation<Rank>
::ApplyConfigSection(const ConfigSection &config) 
{
	//Do this manually for each sub-representation

	if (config.HasValue("innerproduct_algorithm"))
	{
		config.Get("innerproduct_algorithm", Algorithm);
	}
}

template<int Rank> cplx CombinedRepresentation<Rank>
::InnerProduct(const Wavefunction<Rank>& w1, const Wavefunction<Rank>& w2)
{
	blitz::Array<cplx, Rank> d1(w1.GetData());
	blitz::Array<cplx, Rank> d2(w2.GetData());

	/*
	 * Algorithm1 is faster for orthogonal basises
	 * Algorithm2 is faster for nonorthogonal basies
	 *
	 * For a combination of orthogonal and non-orthogonal
	 * basises, they are most likely almost equally fast 
	 */

	if (Algorithm == 1)
	{
		return InnerProductImpl_Algo1(d1, d2);
	}
	else if (Algorithm == 2)
	{
		return InnerProductImpl_Algo2(d1, d2);
	}
	else
	{
		cout << "Unknown InnerProduct algorithm " << Algorithm << endl;
		throw std::runtime_error("Unknown InnerProduct algorithm");
	}

}

/*----------------------------------------------------------------------------
                Implementation of the inner product
  ----------------------------------------------------------------------------*/

template<int Rank> 
cplx CombinedRepresentation<Rank>::InnerProductImpl_Algo1(DataArray d1, DataArray d2)
{
	/* 
	 * Specialized versions of the inner product algorithms are generated
	 * by combinedrepresentation_generator.py and compiled separately
	 */	
	throw std::runtime_error("InnerProduct Not Implemented");
}

template<int Rank> 
cplx CombinedRepresentation<Rank>::InnerProductImpl_Algo2(DataArray d1, DataArray d2)
{
	/* 
	 * Specialized versions of the inner product algorithms are generated
	 * by combinedrepresentation_generator.py and compiled separately
	 */
	throw std::runtime_error("InnerProduct Not Implemented");
}


template class CombinedRepresentation<1>;
template class CombinedRepresentation<2>;
template class CombinedRepresentation<3>;
template class CombinedRepresentation<4>;


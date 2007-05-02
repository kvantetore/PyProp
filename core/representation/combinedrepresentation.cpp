#include "combinedrepresentation.h"

#include "cartesianrepresentation.h"

/*----------------------------------------------------------------------------
                Implementation of CreateSlicedWavefunction
  ----------------------------------------------------------------------------*/

/** General implementation of CreateSlicedWavefunction. This does only throw an error 
  * because I have not been able to write a general algorithm that does this in a
  * Clean way.
  * Insted it should be specialized for every interesting rank
  */
template<int Rank>
Wavefunction<1>* CombinedRepresentation<Rank>
::CreateSlicedWavefunction(const Wavefunction<Rank> &psi, int sliceRank, const blitz::TinyVector<int, Rank-1> &slicePosition)
{
	std::cout << "Invalid rank " << sliceRank << std::endl;
	throw std::runtime_error("Create sliced wavefunction is not implemented for this rank");
}

/** One is not an interesting rank for a combined representation, so the first one  
  * is Rank==2. 
  */
template<> Wavefunction<1>* CombinedRepresentation<2>
::CreateSlicedWavefunction(const Wavefunction<2> &psi, int sliceRank, const blitz::TinyVector<int, 1> &slicePosition)
{
	Wavefunction<1> *slicedPsi = new Wavefunction<1>();
	if (sliceRank == 0)
	{
		slicedPsi->Data.reference(psi.Data(blitz::Range::all(), slicePosition(0)));
	} 
	else
	{
		slicedPsi->Data.reference(psi.Data(slicePosition(0), blitz::Range::all()));
	}
	return slicedPsi;
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



template class CombinedRepresentation<2>;
template class CombinedRepresentation<3>;


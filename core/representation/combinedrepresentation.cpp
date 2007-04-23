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


template class CombinedRepresentation<2>;
template class CombinedRepresentation<3>;


#include "combinedrepresentation.h"
#include "cartesianrepresentation.h"


template<int Rank>
CombinedRepresentation<Rank>::CombinedRepresentation() 
{
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
GetGlobalOverlapMatrix(int rank)
{
	return GetRepresentation(rank)->GetGlobalOverlapMatrix(rank);
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
}

template<int Rank> cplx CombinedRepresentation<Rank>
::InnerProduct(const Wavefunction<Rank>& w1, const Wavefunction<Rank>& w2)
{
	blitz::Array<cplx, Rank> d1(w1.GetData());
	blitz::Array<cplx, Rank> d2(w2.GetData());

	blitz::TinyVector< blitz::Array<double, 1>, Rank> weights;
	blitz::TinyVector< blitz::Array<double, 2>, Rank> overlap;
	blitz::TinyVector< int, Rank> bandWidth;
	bandWidth = 1;
	for (int i=0; i<Rank; i++)
	{
		//Reference all weights
		weights(i).reference(GetLocalWeights(i));

		//Reference all overlap matrices
		if (!this->IsOrthogonalBasis(i))
		{
			bandWidth(i) = GetOverlapBandwidth(i);
			if (this->GetDistributedModel()->IsDistributedRank(i))
			{
				cout << "Rank " << i << " is distributed, and has overlap bandwidth > 1" << endl;
				throw std::runtime_error("Distributed rank has overlap bandwidth > 1");
			}
			
			overlap(i).reference(GetGlobalOverlapMatrix(i));
		}
	}

	double weight = 1;
	cplx innerProduct = 0;
	
	typename blitz::Array<cplx, Rank>::iterator it1 = d1.begin();
	for (int linearCount=0; linearCount<w1.Data.size(); linearCount++)
	{
		weight = 1.;
		for (int curRank=0; curRank<Rank; curRank++)
		{
			weight *= weights(curRank)(it1.position()(curRank));
		}

		/*
		 * Non-orthogonal basises makes inner products more complicated
		 * sum( conj(psi1(i)) * psi2(j) * overlap(i, j) * weight(i), i, j)
		 *
		 * compared to inner product for orthogonal basises
		 * sum( conj(psi1(i)) * psi2(i) * weight(i), i) 
		 *
		 * Overlap is a banded matrix.
		 *
		 * bandWidth == 1 <=> Orthogonal basis
		 *
		 * Note that the overlap matrix is assumed to be stored on BLAS form 
		 * as used by the routine zgbsv.
		 */
		cplx curValue = conj(*it1);

		//Calculate overlap along each rank
		cplx subInnerProduct = 0;
		for (int curRank=0; curRank<Rank; curRank++)
		{
			int curBandWidth = bandWidth(curRank);
			if (curBandWidth > 1)
			{
				blitz::TinyVector<int, Rank> otherIndex = it1.position();
				int curRankIndex = otherIndex(curRank);
				int curRankSize = d1.extent(curRank);
			
				int startBand = std::max(0, curRankIndex - (bandWidth(curRank) - 1) / 2);
				int stopBand = std::min(curRankSize-1, curRankIndex + (bandWidth(curRank) - 1) / 2);
				for (int band=startBand; band<curRankIndex; band++)
				{
					//BLAS index map from "normal" indices
					int Jb = band;
					int Ib = (curRankIndex - band);

					otherIndex(curRank) = band;
					double curOverlap = overlap(curRank)(Jb, Ib); 
					subInnerProduct += curValue * d2(otherIndex) * curOverlap;
				}

				for (int band=curRankIndex; band<=stopBand; band++)
				{
					//BLAS index map from "normal" indices
					int Jb = curRankIndex;
					int Ib = band - curRankIndex;

					otherIndex(curRank) = band;
					double curOverlap = overlap(curRank)(Jb, Ib); 
					subInnerProduct += curValue * d2(otherIndex) * curOverlap;
				}
			}
			else
			{
				subInnerProduct = curValue * d2(it1.position());
			}
		}

		innerProduct += subInnerProduct * weight;
		
		it1++;
	}

	return innerProduct;
}


template class CombinedRepresentation<1>;
template class CombinedRepresentation<2>;
template class CombinedRepresentation<3>;
template class CombinedRepresentation<4>;


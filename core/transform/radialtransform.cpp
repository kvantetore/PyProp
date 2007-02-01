#include "radialtransform.h"
#include "fouriertransform.h"
#include "../representation/representation.h"

using namespace blitz;

template<int Rank>
void RadialTransform<Rank>::TransformRank(Wavefunction<Rank> &psi, int rank, int direction)
{
	if (psi.GetRepresentation().GetDistributedModel().IsDistributedRank(psi, rank))
	{
		std::cout << "Cannot execute fourier transform along distributed rank." << std::endl;
		throw std::runtime_error("Cannot execute fourier transform along distributed rank.");
	}
	
	int minStrideRank = psi.Data.ordering(0);
	int maxStrideRank = psi.Data.ordering(Rank - 1);

	//if we are transforming the minimum strided rank, we can do it better
	//than the generic routine
	if (rank == minStrideRank)
	{
		//std::cout << "Transforming min stride" << std::endl;
		FftOnlyMinStride(psi.Data, direction);
	}
	else if (rank == maxStrideRank)
	{
		//std::cout << "Transforming max stride" << std::endl;
		FftOnlyMaxStride(psi.Data, direction);
	}
	else
	{
		if (rank < Rank/2)
		{
			FftRankNegative(psi.Data, rank, direction);
		} 
		else
		{
			FftRankPositive(psi.Data, rank, direction);
		}
	}	
}

template <int Rank>
void RadialTransform<Rank>::ForwardTransform(Wavefunction<Rank> &psi)
{
	TransformRank(psi, 0, FFT_FORWARD);
}

template <int Rank>
void RadialTransform<Rank>::InverseTransform(Wavefunction<Rank> &psi)
{
	TransformRank(psi, 0, FFT_BACKWARD);
	FftScaleRank(psi, 0);
}


template class RadialTransform<2>;



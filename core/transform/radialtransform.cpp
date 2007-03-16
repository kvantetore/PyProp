#include "radialtransform.h"
#include "fouriertransform.h"
#include "../representation/representation.h"

using namespace blitz;

template<int Rank>
void RadialTransform<Rank>::TransformRank(Wavefunction<Rank> &psi, int rank, int direction)
{
	if (psi.GetRepresentation().GetDistributedModel().IsDistributedRank(rank))
	{
		std::cout << "Cannot execute fourier transform along distributed rank." << std::endl;
		throw std::runtime_error("Cannot execute fourier transform along distributed rank.");
	}
	FftRank(psi.Data, rank, direction);	
}

template <int Rank>
void RadialTransform<Rank>::ForwardTransform(Wavefunction<Rank> &psi, int rank)
{
	TransformRank(psi, rank, FFT_FORWARD);
}

template <int Rank>
void RadialTransform<Rank>::InverseTransform(Wavefunction<Rank> &psi, int rank)
{
	TransformRank(psi, rank, FFT_BACKWARD);
	FftScaleRank(psi.Data, rank);
}


template class RadialTransform<2>;
template class RadialTransform<3>;



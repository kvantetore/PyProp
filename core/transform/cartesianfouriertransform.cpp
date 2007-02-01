#include "cartesianfouriertransform.h"
#include "../representation/representation.h"
#include "../representation/cartesianrepresentation.h"

template<int Rank>
void CartesianFourierTransform<Rank>::TransformExceptDistributedRank(Wavefunction<Rank> &psi, int direction)
{
	if (!psi.GetRepresentation().GetDistributedModel().HasDistributedRangeMaxStride(psi))
	{
		throw std::runtime_error("Maximum stride is not distributed, this case is not implemented in cartesian fourier transform.");
	}
	FftAllExceptMaxStride(psi.Data, FFT_FORWARD);
}

template<int Rank>
void CartesianFourierTransform<Rank>::TransformRank(Wavefunction<Rank> &psi, int rank, int direction)
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

template<>
void CartesianFourierTransform<1>::TransformRank(Wavefunction<1> &psi, int rank, int direction)
{
	FftRank(psi.Data, rank, direction);
}


template<int Rank>
void CartesianFourierTransform<Rank>::Renormalize(Wavefunction<Rank> &psi)
{
	FftScale(psi);
}

template<int Rank>
void CartesianFourierTransform<Rank>::ForwardTransform(Wavefunction<Rank> &psi)
{
	FourierTransform(psi, FFT_FORWARD);
}

template<int Rank>
void CartesianFourierTransform<Rank>::InverseTransform(Wavefunction<Rank> &psi)
{
	FourierTransform(psi, FFT_BACKWARD);
}

/** Create representations for a full transform */
template<int Rank> 
CartesianRepresentation<Rank> CartesianFourierTransform<Rank>::CreateFourierRepresentation(const CartesianRepresentation<Rank> &gridRepresentation)
{
	//A fourier representation is basically another cartesian representation 
	//(of the same size) with a different grid
	CartesianRepresentation<Rank> fftRepr(gridRepresentation);
	
	for (int i=0;i<Rank;i++)
	{
		double kMin = - M_PI / gridRepresentation.Range(i).Dx;
		double kMax = - kMin;
		//std::cout << " [" << kMin << ", " << kMax << ", " << Range(i).Count << "]" << std::endl;
		fftRepr.Range(i) = CartesianRange(kMin, kMax, gridRepresentation.Range(i).Count, true);
	}
	return fftRepr;
}		
	
/** Create representations for partial transforms */
template<int Rank>
CartesianRepresentation<Rank> CartesianFourierTransform<Rank>::CreateFourierRepresentation(const CartesianRepresentation<Rank> &gridRepresentation, int rank)
{
	CartesianRepresentation<Rank> fftRepr(gridRepresentation);
	
	double kMin = - M_PI / gridRepresentation.Range(rank).Dx;
	double kMax = - kMin;
	fftRepr.Range(rank) = CartesianRange(kMin, kMax, gridRepresentation.Range(rank).Count, true);
	
	return fftRepr;
}


template class CartesianFourierTransform<1>;
template class CartesianFourierTransform<2>;
template class CartesianFourierTransform<3>;
template class CartesianFourierTransform<4>;
template class CartesianFourierTransform<5>;
template class CartesianFourierTransform<6>;
template class CartesianFourierTransform<7>;
template class CartesianFourierTransform<8>;
template class CartesianFourierTransform<9>;

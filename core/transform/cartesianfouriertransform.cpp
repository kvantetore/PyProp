#include "cartesianfouriertransform.h"
#include "../representation/representation.h"
#include "../representation/cartesianrepresentation.h"

template<int Rank>
void CartesianFourierTransform<Rank>::TransformRank(Wavefunction<Rank> &psi, int rank, int direction)
{
	if (psi.GetRepresentation()->GetDistributedModel()->IsDistributedRank(rank))
	{
		std::cout << "Cannot execute fourier transform along distributed rank." << std::endl;
		throw std::runtime_error("Cannot execute fourier transform along distributed rank.");
	}

	FftRank(psi.Data, rank, direction);
}

/* Specialized TransformRank for 1D */
template<>
void CartesianFourierTransform<1>::TransformRank(Wavefunction<1> &psi, int rank, int direction)
{
	FftRank(psi.Data, rank, direction);
}

/*
 * After a performing forward and backward transforms with one rank at a time,
 * a call to Renormalize() is required to rescale the wavefunction to the expected
 * norm. This is done automatically by ForwardTransform and InverseTransform
 */
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


template<int Rank>
void CartesianFourierTransform<Rank>::FourierTransform(Wavefunction<Rank> &psi, int direction)
{
	if (psi.GetRepresentation()->GetDistributedModel()->ProcCount == 1)
	{
		//For one processor, we can execute a full Rank-D FFT 
		FftAll(psi.Data, direction);
	}
	else 
	{
		throw std::runtime_error("Cannot execute fourier transform along distributed rank.");
	}
	
	if (direction == FFT_BACKWARD) 
	{
		FftScale(psi);		
	}
}

/** Create representations for a full transform */
template<int Rank> 
typename CartesianRepresentation<Rank>::Ptr CartesianFourierTransform<Rank>::CreateFourierRepresentation(const CartesianRepresentation<Rank> &gridRepresentation)
{
	//A fourier representation is basically another cartesian representation 
	//(of the same size) with a different grid
	typedef CartesianRepresentation<Rank> Repr;
	typedef typename CartesianRepresentation<Rank>::Ptr ReprPtr;
	ReprPtr fftRepr(new Repr(gridRepresentation));

	for (int i=0;i<Rank;i++)
	{
		double kMin = - M_PI / gridRepresentation.Range(i).Dx;
		double kMax = - kMin;
		//std::cout << " [" << kMin << ", " << kMax << ", " << Range(i).Count << "]" << std::endl;
		fftRepr->Range(i) = CartesianRange(kMin, kMax, gridRepresentation.Range(i).Count, true);
	}
	return fftRepr;
}		
	
/** Create representations for partial transforms */
template<int Rank>
typename CartesianRepresentation<Rank>::Ptr CartesianFourierTransform<Rank>::CreateFourierRepresentation(const CartesianRepresentation<Rank> &gridRepresentation, int rank)
{
	typedef CartesianRepresentation<Rank> Repr;
	typedef typename CartesianRepresentation<Rank>::Ptr ReprPtr;
	ReprPtr fftRepr(new Repr(gridRepresentation));
	
	double kMin = - M_PI / gridRepresentation.Range(rank).Dx;
	double kMax = - kMin;
	fftRepr->Range(rank) = CartesianRange(kMin, kMax, gridRepresentation.Range(rank).Count, true);
	
	return fftRepr;
}

/**
Scales the wavefunction after a full forward / backward transform.
It is equivalent to calling FftScaleRank() for all ranks (only much
more efficient).
*/
template <int Rank>
void FftScale(Wavefunction<Rank> &psi)
{
	double scale = 1.0;
	blitz::TinyVector<int, Rank> fullShape = psi.GetRepresentation()->GetFullShape();
	for (int dimension=0; dimension<Rank; dimension++)
	{
		scale /= fullShape(dimension);
	}
	psi.Data = psi.Data * scale;		
}

/**
Scales the wavefunction after a forward / backward transform along
the specified rank. 
*/
template <int Rank>
void FftScaleRank(Wavefunction<Rank> &psi, int rank)
{
	double scale = 1.0;
	blitz::TinyVector<int, Rank> fullShape = psi.GetRepresentation()->GetFullShape();
	scale /= fullShape(rank);
	psi.Data = psi.Data * scale;		
}


template class CartesianFourierTransform<1>;
template class CartesianFourierTransform<2>;
template class CartesianFourierTransform<3>;
template class CartesianFourierTransform<4>;





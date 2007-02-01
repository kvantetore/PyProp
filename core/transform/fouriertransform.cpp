#include "fouriertransform.h"

#include <complex>
#include <fftw3.h>

#include "../common.h"
#include "../index.h"
#include "../wavefunction.h"
#include "../representation/representation.h"
#include "../mpi/distributedmodel.h"

#define FFT_FLAG FFTW_MEASURE

/* ----------------------------------------------------------------------- */
/* ----------              Low level routines                   ---------- */
/* ----------------------------------------------------------------------- */

/**
Transforms all ranks except the one with the maximum stride
(usually conforming to the distributed rank)
*/
template<int Rank>
void FftAllExceptMaxStride(blitz::Array<cplx, Rank> &array, int direction)
{
	fftw_complex *data = reinterpret_cast<fftw_complex*>( array.data() );
	
	blitz::TinyVector<int, Rank-1> shape;
	int fftDistance = 1;
	for (int i=0; i<Rank-1; i++)
	{
		shape(i) = array.extent(array.ordering(Rank - i - 2));
		fftDistance *= shape(i);
	}
	int fftCount = array.extent(array.ordering(Rank-1));
	int fftStride = 1;
	
	/*
	std::cout << Rank-1 << "D FFT, count = " << fftCount 
	          << ", distance = " << fftDistance 
		  << std::endl; 
	*/

	fftw_plan plan = fftw_plan_many_dft(
				Rank-1,		//Rank
				&shape(0), 	//N[Rank]
				fftCount,	//How many
				data, 		//In-ptr
				0, 		//In-LogicalSize[Rank] (IGNORED)
				fftStride, 	//In-stride
				fftDistance,	//In-Distance
				data, 		//Out-ptr
				0,		//Out-LogicalSze[Rank] (IGNORED)
				fftStride, 	//Out-stride
				fftDistance,	//Out-distance
				direction, 	//Transform direction
				FFT_FLAG	//What kind of a plan
			);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}

/* 
Specialication of FftAllExceptMaxStride for Rank == 1
*/
template<>
void FftAllExceptMaxStride(blitz::Array<cplx, 1> &array, int direction)
{
	//do nothing for Rank == 1
}


/**
Transforms all ranks in an array. This only works when
none of the ranks are distributed
*/
template<int Rank>
void FftAll(blitz::Array<cplx, Rank> &array, int direction)
{
	fftw_complex *data = reinterpret_cast<fftw_complex*>( array.data() );
	
	blitz::TinyVector<int, Rank> shape;
	for (int i=0; i<Rank; i++)
	{
		shape(i) = array.extent(array.ordering(Rank - i - 1));
		//std::cout << "stride(" << i << ") = " << array.stride(array.ordering(Rank - i - 1));
	}
	
	/*
	std::cout << Rank-1 << "D FFT, count = " << fftCount 
	          << ", distance = " << fftDistance 
		  << std::endl; 
	*/

	fftw_plan plan = fftw_plan_many_dft(
				Rank,		//Rank
				&shape(0), 	//N[Rank]
				1,		//How many
				data, 		//In-ptr
				0, 		//In-LogicalSize[Rank] (IGNORED)
				1, 		//In-stride
				1,		//In-Distance
				data, 		//Out-ptr
				0,		//Out-LogicalSze[Rank] (IGNORED)
				1,	 	//Out-stride
				1,		//Out-distance
				direction, 	//Transform direction
				FFT_FLAG	//What kind of a plan
			);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}

/**
Transforms the direction with the minimum stride.
After a change of distribution, this will usually be the
rank that was previously distributed.
*/
template<int Rank>
void FftOnlyMinStride(blitz::Array<cplx, Rank> &array, int direction)
{
	fftw_complex *data = reinterpret_cast<fftw_complex*>( array.data() );
	
	int fftDistance = array.extent(array.ordering(0));;
	int fftStride = 1;
	int fftSize = array.extent(array.ordering(0));;
	int fftCount = array.size() / fftSize;

	fftw_plan plan = fftw_plan_many_dft(
				1,		//Rank
				&fftSize, 	//N[Rank]
				fftCount,	//How many
				data, 		//In-ptr
				0, 		//In-LogicalSize[Rank] (IGNORED)
				fftStride, 	//In-stride
				fftDistance,	//In-Distance
				data, 		//Out-ptr
				0,		//Out-LogicalSze[Rank] (IGNORED)
				fftStride, 	//Out-stride
				fftDistance,	//Out-distance
				direction, 	//Transform direction
				FFT_FLAG	//What kind of a plan
			);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}


/**
Transforms the direction with the minimum stride.
After a change of distribution, this will usually be the
rank that was previously distributed.
*/
template<int Rank>
void FftOnlyMaxStride(blitz::Array<cplx, Rank> &array, int direction)
{
	fftw_complex *data = reinterpret_cast<fftw_complex*>( array.data() );
	
	int fftDistance = 1;
	int fftStride = array.stride(array.ordering(Rank - 1));
	int fftSize = array.extent(array.ordering(Rank - 1));
	int fftCount = array.size() / fftSize;

	fftw_plan plan = fftw_plan_many_dft(
				1,		//Rank
				&fftSize, 	//N[Rank]
				fftCount,	//How many
				data, 		//In-ptr
				0, 		//In-LogicalSize[Rank] (IGNORED)
				fftStride, 	//In-stride
				fftDistance,	//In-Distance
				data, 		//Out-ptr
				0,		//Out-LogicalSze[Rank] (IGNORED)
				fftStride, 	//Out-stride
				fftDistance,	//Out-distance
				direction, 	//Transform direction
				FFT_FLAG	//What kind of a plan
			);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}


/**
Transforms the array only along the specified rank using the "Positive" algorithm 
(generalized algorithm based on FftOnlyMinStride)
*/
template<int Rank>
void FftRankPositive(blitz::Array<cplx, Rank> &array, int rank, int direction)
{
	//Find the order index of our rank
	int rankOrderIndex = -1;
	for (int i=0; i<Rank; i++)
	{
		if (array.ordering(i) == rank)
		{
			rankOrderIndex = i;
			break;
		}
	}
	if (rankOrderIndex < 0)
	{
		std::cout << "WHAT! Could not find rank " << rank 
		          << " in ordering of array of rank " << Rank << std::endl;
		throw std::runtime_error("Error in FftRankPositive");
	}

	//FFT Parameters
	int fftStride = array.stride(rank);
	int fftSize = array.extent(rank);
	int fftDist = 1;
	for (int i=0; i<=rankOrderIndex; i++)
	{
		fftDist *= array.extent(array.ordering(i));
	}
	int fftHowMany = 1;
	for (int i=rankOrderIndex+1; i<Rank; i++)
	{
		fftHowMany *= array.extent(array.ordering(i));
	}
	
	//Repeat parameters
	int repeatDist = 1;
	int repeatCount = 1;
	for (int i=0; i<rankOrderIndex; i++)
	{
		repeatCount *= array.extent(array.ordering(i));
	}

	if (repeatCount > fftHowMany)
	{
		std::cout << "Warning: Using FftRankPositive when it would be more "
			  << "efficient to use FftRankNegative." << std::endl
			  << "Consider using FftRankNegative instead. " << std::endl
			  << "  repeatCount = " << repeatCount << std::endl
			  << "  fftCount    = " << fftHowMany << std::endl;			  
	}
		
	fftw_complex *data = reinterpret_cast<fftw_complex*>( array.data() );
	//Create a plan for the transform of ONE array
	fftw_plan plan = fftw_plan_many_dft(
				1, 		//Rank
				&fftSize, 	//N[Rank]
				fftHowMany,	//How many
				data, 		//In-ptr
				0, 		//In-LogicalSize[Rank] (IGNORED)
				fftStride, 	//In-stride
				fftDist,	//In-Distance (IGNORED)
				data, 		//Out-ptr
				0,		//Out-LogicalSze[Rank] (IGNORED)
				fftStride, 	//Out-stride
				fftDist,	//Out-distance (IGNORED)
				direction, 	//Transform direction
				FFT_FLAG	//What kind of a plan
			);
	
	
	for (int repeatIndex=0; repeatIndex<repeatCount; repeatIndex++)
	{
		fftw_execute_dft(plan, data, data);
		data += repeatDist;
	}
	fftw_destroy_plan(plan);
	
}


/**
Transforms the array only along the specified rank using the "Negative" algorithm 
(generalized algorithm based on FftOnlyMaxStride)
*/
template<int Rank>
void FftRankNegative(blitz::Array<cplx, Rank> &array, int rank, int direction)
{
	//Find the order index of our rank
	int rankOrderIndex = -1;
	for (int i=0; i<Rank; i++)
	{
		if (array.ordering(i) == rank)
		{
			rankOrderIndex = i;
			break;
		}
	}
	if (rankOrderIndex < 0)
	{
		std::cout << "WHAT! Could not find rank " << rank 
		          << " in ordering of array of rank " << Rank << std::endl;
		throw std::runtime_error("Error in FftRankPositive");
	}

	//FFT Parameters
	int fftStride = array.stride(rank);
	int fftSize = array.extent(rank);
	int fftDist = 1;
	int fftHowMany = 1;
	for (int i=0; i<rankOrderIndex; i++)
	{
		fftHowMany *= array.extent(array.ordering(i));
	}
	
	//Repeat parameters
	int repeatDist = 1;
	for (int i=0; i<=rankOrderIndex; i++)
	{
		repeatDist *= array.extent(array.ordering(i));
	}
	int repeatCount = 1;
	for (int i=rankOrderIndex+1; i<Rank; i++)
	{
		repeatCount *= array.extent(array.ordering(i));
	}
	
	if (repeatCount > fftHowMany)
	{
		std::cout << "Warning: Using FftRankNegative when it would be more "
			  << "efficient to use FftRankPositive." << std::endl
			  << "Consider using FftRankPositive instead. " << std::endl
			  << "  repeatCount = " << repeatCount << std::endl
			  << "  fftCount    = " << fftHowMany << std::endl;			  
	}
	
	fftw_complex *data = reinterpret_cast<fftw_complex*>( array.data() );
	//Create a plan for the transform of ONE array
	fftw_plan plan = fftw_plan_many_dft(
				1, 		//Rank
				&fftSize, 	//N[Rank]
				fftHowMany,	//How many
				data, 		//In-ptr
				0, 		//In-LogicalSize[Rank] (IGNORED)
				fftStride, 	//In-stride
				fftDist,	//In-Distance (IGNORED)
				data, 		//Out-ptr
				0,		//Out-LogicalSze[Rank] (IGNORED)
				fftStride, 	//Out-stride
				fftDist,	//Out-distance (IGNORED)
				direction, 	//Transform direction
				FFT_FLAG	//What kind of a plan
			);
	
	for (int repeatIndex=0; repeatIndex<repeatCount; repeatIndex++)
	{
		fftw_execute_dft(plan, data, data);
		data += repeatDist;
	}
	fftw_destroy_plan(plan);
}


/**
Transforms the array only along the specified rank
*/
template<int Rank>
void FftRank(blitz::Array<cplx, Rank> &array, int rank, int direction)
{
	int count = array.extent(rank);
	int stride = array.stride(rank);
	
	blitz::TinyVector<int, Rank> index = 0;
	blitz::TinyVector<int, Rank> extent = array.shape();
	
	do
	{
		fftw_complex *data = reinterpret_cast<fftw_complex*>( & array(index) );
	
		//Create a plan for the transform of ONE array
		fftw_plan plan = fftw_plan_many_dft(
					1, 		//Rank
					&count, 	//N[Rank]
					1,		//How many
					data, 		//In-ptr
					0, 		//In-LogicalSize[Rank] (IGNORED)
					stride, 	//In-stride
					1,		//In-Distance (IGNORED)
					data, 		//Out-ptr
					0,		//Out-LogicalSze[Rank] (IGNORED)
					stride, 	//Out-stride
					1,		//Out-distance (IGNORED)
					direction, 	//Transform direction
					FFT_FLAG	//What kind of a plan
				);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
	
	} while (IncIndex(extent, index, Rank - 1, rank));
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
	blitz::TinyVector<int, Rank> fullShape = psi.GetRepresentation().GetFullShape();
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
	blitz::TinyVector<int, Rank> fullShape = psi.GetRepresentation().GetFullShape();
	scale /= fullShape(rank);
	psi.Data = psi.Data * scale;		
}

/* ----------------------------------------------------------------------- */
/* ----------              High level routines                  ---------- */
/* ----------------------------------------------------------------------- */
template<int Rank>
void FourierTransform(Wavefunction<Rank> &psi, int direction)
{
	if (psi.GetRepresentation().GetDistributedModel().ProcCount == 1)
	{
		//For one processor, we can execute a full Rank-D FFT 
		FftAll(psi.Data, direction);
	}
	else 
	{
		//On a multi processor system, first transform all but the distributed rank,
		//change distribution, and hopefully, the previously distrubted rank is now
		//the min stride rank.
		if (!psi.GetRepresentation().GetDistributedModel().HasDistributedRangeMaxStride(psi))
		{
			throw std::runtime_error("Error in FourierTransform(psi, direction): Maximum stride is not distributed, this case is not implemented in FourierTransform() yet.");
		}
		int distributedRank = psi.GetRepresentation().GetDistributedModel().GetDistributedRank(psi);
		FftAllExceptMaxStride(psi.Data, direction);
		
		//change representation
		psi.GetRepresentation().GetDistributedModel().ChangeRepresentation(psi);
		
		int minStrideRank = psi.Data.ordering(0);
		//if we are transforming the minimum strided rank, we can do it better
		//than the generic routine
		if (distributedRank == minStrideRank)
		{
			FftOnlyMinStride(psi.Data, direction);
		}
		else
		{
			std::cout << "Warning from FourierTransform(psi,direction): Performing FFT along a suboptimal rank " << distributedRank << std::endl;
			FftRank(psi.Data, distributedRank, direction);
		}
	}
	
	if (direction == FFTW_BACKWARD) 
	{
		FftScale(psi);		
	}
}

template<int Rank>
void FourierTransformDimension(Wavefunction<Rank> &psi, int dimension, int direction)
{
	if (psi.GetRepresentation().GetDistributedModel().GetDistributedRank(psi) == dimension)
	{
		std::cout << "Rank " << dimension << " is currently distributed, and fft cannot be computed "
		          << "be computed on a distributed rank" << std::endl;
		throw std::runtime_error("Error in fourier transform");
	}
	
	FftRank(psi.Data, dimension, direction);
}


template void FftAll(blitz::Array<cplx, 1> &array, int direction);
//explicitly specialized:
//template void FftAllExceptMaxStride(blitz::Array<cplx, 1> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 1> &array, int direction);
template void FftAll(blitz::Array<cplx, 2> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 2> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 2> &array, int direction);
template void FftAll(blitz::Array<cplx, 3> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 3> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 3> &array, int direction);
template void FftAll(blitz::Array<cplx, 4> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 4> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 4> &array, int direction);
template void FftAll(blitz::Array<cplx, 5> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 5> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 5> &array, int direction);
template void FftAll(blitz::Array<cplx, 6> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 6> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 6> &array, int direction);
template void FftAll(blitz::Array<cplx, 7> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 7> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 7> &array, int direction);
template void FftAll(blitz::Array<cplx, 8> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 8> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 8> &array, int direction);
template void FftAll(blitz::Array<cplx, 9> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 9> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 9> &array, int direction);

template void FftOnlyMaxStride(blitz::Array<cplx, 1> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 2> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 3> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 4> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 5> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 6> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 7> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 8> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 9> &array, int direction);

template void FftRank(blitz::Array<cplx, 1> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 2> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 3> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 4> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 5> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 6> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 7> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 8> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 9> &array, int rank, int direction);

template void FftRankPositive(blitz::Array<cplx, 1> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 2> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 3> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 4> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 5> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 6> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 7> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 8> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 9> &array, int rank, int direction);

template void FftRankNegative(blitz::Array<cplx, 1> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 2> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 3> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 4> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 5> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 6> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 7> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 8> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 9> &array, int rank, int direction);

template void FftScale(Wavefunction<1> &psi);
template void FftScale(Wavefunction<2> &psi);
template void FftScale(Wavefunction<3> &psi);
template void FftScale(Wavefunction<4> &psi);
template void FftScale(Wavefunction<5> &psi);
template void FftScale(Wavefunction<6> &psi);
template void FftScale(Wavefunction<7> &psi);
template void FftScale(Wavefunction<8> &psi);
template void FftScale(Wavefunction<9> &psi);

template void FftScaleRank(Wavefunction<1> &psi, int rank);
template void FftScaleRank(Wavefunction<2> &psi, int rank);
template void FftScaleRank(Wavefunction<3> &psi, int rank);
template void FftScaleRank(Wavefunction<4> &psi, int rank);
template void FftScaleRank(Wavefunction<5> &psi, int rank);
template void FftScaleRank(Wavefunction<6> &psi, int rank);
template void FftScaleRank(Wavefunction<7> &psi, int rank);
template void FftScaleRank(Wavefunction<8> &psi, int rank);
template void FftScaleRank(Wavefunction<9> &psi, int rank);

template void FourierTransform(Wavefunction<1> &psi, int direction);
template void FourierTransform(Wavefunction<2> &psi, int direction);
template void FourierTransform(Wavefunction<3> &psi, int direction);
template void FourierTransform(Wavefunction<4> &psi, int direction);
template void FourierTransform(Wavefunction<5> &psi, int direction);
template void FourierTransform(Wavefunction<6> &psi, int direction);
template void FourierTransform(Wavefunction<7> &psi, int direction);
template void FourierTransform(Wavefunction<8> &psi, int direction);
template void FourierTransform(Wavefunction<9> &psi, int direction);

template void FourierTransformDimension(Wavefunction<1> &psi, int dimension, int direction);
template void FourierTransformDimension(Wavefunction<2> &psi, int dimension, int direction);
template void FourierTransformDimension(Wavefunction<3> &psi, int dimension, int direction);
template void FourierTransformDimension(Wavefunction<4> &psi, int dimension, int direction);
template void FourierTransformDimension(Wavefunction<5> &psi, int dimension, int direction);
template void FourierTransformDimension(Wavefunction<6> &psi, int dimension, int direction);
template void FourierTransformDimension(Wavefunction<7> &psi, int dimension, int direction);
template void FourierTransformDimension(Wavefunction<8> &psi, int dimension, int direction);
template void FourierTransformDimension(Wavefunction<9> &psi, int dimension, int direction);

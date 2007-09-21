#include "fouriertransform.h"

#include <complex>
#include <fftw3.h>

#define FFT_FLAG FFTW_ESTIMATE
//#define FFT_FLAG FFTW_MEASURE <--- Dangerous.


/* --------------- General Routines -----------------------  */
/*
 * These routines work on one rank at a time, but is much more efficient than
 * looping over all other ranks and perform 1D fft (because fftw is more optimized
 * than general c++ code)
 *
 */


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
	int minStrideRank = array.ordering(0);
	int maxStrideRank = array.ordering(Rank - 1);

	//if we are transforming the minimum strided rank, we can do it better
	//than the generic routine
	if (rank == minStrideRank)
	{
		FftOnlyMinStride(array, direction);
	}
	else if (rank == maxStrideRank)
	{
		FftOnlyMaxStride(array, direction);
	}
	else
	{
		int rankOrderIndex = -1;
		for (int i=0; i<Rank; i++)
		{
			if (array.ordering(i) == rank)
			{
				rankOrderIndex = i;
				break;
			}
		}
	
		int lowerCount = 1;
		for (int i=0; i<rankOrderIndex; i++)
		{
			lowerCount *= array.extent(array.ordering(i));
		}
		int upperCount = 1;
		for (int i=rankOrderIndex+1; i<Rank; i++)
		{
			upperCount *= array.extent(array.ordering(i));
		}
		if (lowerCount > upperCount)
		{
			FftRankNegative(array, rank, direction);
		} 
		else
		{
			FftRankPositive(array, rank, direction);
		}
	}	
}

/* --------------- Specialized functions --------------- */
/* These functions provide a more specialized implementation of fourier transform
 * along all, max stride, or min stride rank. They are more efficient than transforming
 * one rank at a time.
 */

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
				1,			//Rank
				&fftSize, 	//N[Rank]
				fftCount,	//How many
				data, 		//In-ptr
				0, 			//In-LogicalSize[Rank] (IGNORED)
				fftStride, 	//In-stride
				fftDistance,//In-Distance
				data, 		//Out-ptr
				0,			//Out-LogicalSze[Rank] (IGNORED)
				fftStride, 	//Out-stride
				fftDistance,//Out-distance
				direction, 	//Transform direction
				FFT_FLAG	//What kind of a plan
			);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}

template<int Rank>
void FftScaleRank(blitz::Array<cplx, Rank> &array, int rank)
{
	double scale = array.extent(rank);
	array = array / scale;
}


template void FftAll(blitz::Array<cplx, 1> &array, int direction);
template void FftAll(blitz::Array<cplx, 2> &array, int direction);
template void FftAll(blitz::Array<cplx, 3> &array, int direction);
template void FftAll(blitz::Array<cplx, 4> &array, int direction);


//explicitly specialized:
//template void FftAllExceptMaxStride(blitz::Array<cplx, 1> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 2> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 3> &array, int direction);
template void FftAllExceptMaxStride(blitz::Array<cplx, 4> &array, int direction);

template void FftOnlyMinStride(blitz::Array<cplx, 1> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 2> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 3> &array, int direction);
template void FftOnlyMinStride(blitz::Array<cplx, 4> &array, int direction);


template void FftOnlyMaxStride(blitz::Array<cplx, 1> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 2> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 3> &array, int direction);
template void FftOnlyMaxStride(blitz::Array<cplx, 4> &array, int direction);


template void FftRank(blitz::Array<cplx, 1> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 2> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 3> &array, int rank, int direction);
template void FftRank(blitz::Array<cplx, 4> &array, int rank, int direction);


template void FftRankPositive(blitz::Array<cplx, 1> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 2> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 3> &array, int rank, int direction);
template void FftRankPositive(blitz::Array<cplx, 4> &array, int rank, int direction);


template void FftRankNegative(blitz::Array<cplx, 1> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 2> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 3> &array, int rank, int direction);
template void FftRankNegative(blitz::Array<cplx, 4> &array, int rank, int direction);



template void FftScaleRank(blitz::Array<cplx, 1> &array, int rank);
template void FftScaleRank(blitz::Array<cplx, 2> &array, int rank);
template void FftScaleRank(blitz::Array<cplx, 3> &array, int rank);
template void FftScaleRank(blitz::Array<cplx, 4> &array, int rank);



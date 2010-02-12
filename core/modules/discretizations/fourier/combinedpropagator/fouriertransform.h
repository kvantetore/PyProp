#ifndef FOURIERTRANSFORM_H
#define FOURIERTRANSFORM_H

#include "../common.h"

#define FFT_FORWARD -1
#define FFT_BACKWARD 1

/*
Low level routines
These routines does _NOT_ have awareness of distributed rank, and it is up to the caller 
to ensure that they don't transform along a distributed rank
*/

//General routines
template<int Rank> void FftRank(blitz::Array<cplx, Rank> &array, int rank, int direction);
template<int Rank> void FftRankPositive(blitz::Array<cplx, Rank> &array, int rank, int direction);
template<int Rank> void FftRankNegative(blitz::Array<cplx, Rank> &array, int rank, int direction);

//Specialized routines
template<int Rank> void FftAll(blitz::Array<cplx, Rank> &array, int direction);
template<int Rank> void FftAllExceptMaxStride(blitz::Array<cplx, Rank> &array, int direction);
template<int Rank> void FftOnlyMinStride(blitz::Array<cplx, Rank> &array, int direction);
template<int Rank> void FftOnlyMaxStride(blitz::Array<cplx, Rank> &array, int direction);

//Scaling
template<int Rank> void FftScaleRank(blitz::Array<cplx, Rank> &array, int rank);

#endif

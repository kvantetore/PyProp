#ifndef FOURIERTRANSFORM_H
#define FOURIERTRANSFORM_H

#include "../common.h"
#include "../wavefunction.h"

#define FFT_FORWARD -1
#define FFT_BACKWARD 1

/*
High level routines
These routines have awareness of distributed rank, and will not attempt to transform
along distributed rank.
They also scale the wavefunction if the direction is FFT_BACKWARD, such that
*/
template<int Rank> void FourierTransformDimension(Wavefunction<Rank> &psi, int dimension, int direction);
template<int Rank> void FourierTransform(Wavefunction<Rank> &psi, int direction);

/*
Low level routines
These routines does _NOT_ have awareness of distributed rank, and it is up to the caller 
to ensure that they don't transform along a distributed rank
*/
template<int Rank> void FftAll(blitz::Array<cplx, Rank> &array, int direction);
template<int Rank> void FftAllExceptMaxStride(blitz::Array<cplx, Rank> &array, int direction);
template<int Rank> void FftOnlyMinStride(blitz::Array<cplx, Rank> &array, int direction);
template<int Rank> void FftOnlyMaxStride(blitz::Array<cplx, Rank> &array, int direction);
template<int Rank> void FftRank(blitz::Array<cplx, Rank> &array, int rank, int direction);
template<int Rank> void FftRankPositive(blitz::Array<cplx, Rank> &array, int rank, int direction);
template<int Rank> void FftRankNegative(blitz::Array<cplx, Rank> &array, int rank, int direction);

template<int Rank> void FftScale(Wavefunction<Rank> &psi);
template<int Rank> void FftScaleRank(Wavefunction<Rank> &psi, int rank);

#endif

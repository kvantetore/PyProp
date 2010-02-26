#ifndef MATRIX_CONVERSION_H
#define MATRIX_CONVERSION_H

#include "../common.h"

blitz::Array<cplx, 2> ConvertMatrixBlasBandedToFull(blitz::Array<cplx, 2> blasBanded);
blitz::Array<cplx, 2> ConvertMatrixBlasBandedToDistributedBanded(blitz::Array<cplx, 2> blasBanded, int localGridSize, int localGridStart);
blitz::Array<cplx, 2> ConvertMatrixBlasBandedToLapackBanded(blitz::Array<cplx, 2> blasBanded);
blitz::Array<cplx, 1> GetDiagonalViewLapackBanded(blitz::Array<cplx, 2> lapackBanded);

#endif


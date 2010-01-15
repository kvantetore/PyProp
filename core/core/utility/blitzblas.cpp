#include "blitzblas.h"

using namespace blitz;

/*
 * Decide which version of blitzblas we should use, the reference
 * implementation, using blitz array expressions, or the much faster
 * blas implementation
 */
#ifndef PYPROP_USE_BLAS
#include "blitzblas_ref.cpp"
#else 
#ifdef PYPROP_USE_BLAS_ACML
#include "blitzblas_acml.cpp"
#else
#include "blitzblas_cblas.cpp"
#endif
#endif  //PYPROP_USE_BLAS


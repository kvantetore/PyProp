#include "blitzblas.h"

using namespace blitz;

/*
 * Decide which version of blitzblas we should use, the reference
 * implementation, using blitz array expressions, or the much faster
 * blas implementation
 */
#ifdef PYPROP_USE_BLAS_ACML
#include "blitzblas_acml.cpp"
#else
#include "blitzblas_cblas.cpp"
#endif

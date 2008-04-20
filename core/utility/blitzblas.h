#ifndef BLITZBLAS_H
#define BLITZBLAS_H

#include "../common.h" 

//Performs the matrix-matrix product C = A B
void MatrixMatrixMultiply(const blitz::Array<cplx, 2> &A, const blitz::Array<cplx, 2> &B, blitz::Array<cplx, 2> &C);
void MatrixMatrixMultiply(const blitz::Array<double, 2> &A, const blitz::Array<double, 2> &B, blitz::Array<double, 2> &C);

//Performs the matrix-vector product w = A v
void MatrixVectorMultiply(const blitz::Array<cplx, 2> &A, const blitz::Array<cplx, 1> &v, blitz::Array<cplx, 1> &w);
void MatrixVectorMultiply(const blitz::Array<double, 2> &A, const blitz::Array<double, 1> &v, blitz::Array<double, 1> &w);

//Performs the matrix-vector product x = A y for hermitian banded A
#ifndef PYPROP_USE_BLAS_ACML
void MatrixVectorMultiplyHermitianBanded(const blitz::Array<cplx, 2> &A, const blitz::Array<cplx, 1> &x, 
	blitz::Array<cplx, 1> &y, cplx alpha, cplx beta);

//Performs the matrix-vector product x = A y for banded A
void MatrixVectorMultiplyBanded(const blitz::Array<cplx, 2> &A, const blitz::Array<cplx, 1> &x, 
	blitz::Array<cplx, 1> &y, cplx alpha, cplx beta, int M);
#endif	

//Performs the vector-vector elementwisse product w_i = u_i * v_i
//any of u, v, and w may be the same vector
template<int Rank>
void VectorElementMultiply(const blitz::Array<cplx, Rank> &u, const blitz::Array<cplx, Rank> &v, blitz::Array<cplx, Rank> &w);
template<int Rank>
void VectorElementMultiply(const blitz::Array<double, Rank> &u, const blitz::Array<double, Rank> &v, blitz::Array<double, Rank> &w);


//Performs the inner product conj(u) * v, and returns the value
//u and v may be the same vector
template<int Rank>
double VectorInnerProduct(const blitz::Array<double, Rank> &u, const blitz::Array<double, Rank> &v);
template<int Rank>
cplx VectorInnerProduct(const blitz::Array<cplx, Rank> &u, const blitz::Array<cplx, Rank> &v);

template<int Rank>
void ScaleVector(cplx scaling, blitz::Array<cplx, Rank> &x)
{
	x *= scaling;
}

template<int Rank>
void CopyVector(cplx sourceScaling, const blitz::Array<cplx, Rank> &source, cplx destScaling, blitz::Array<cplx, Rank> &dest)
{
	dest = sourceScaling * source + destScaling * dest;
}

#endif


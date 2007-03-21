#ifndef BLITZBLAS_H
#define BLITZBLAS_H

#include "../common.h" 

//Performs the matrix-matrix product C = A B
void MatrixMatrixMultiply(const blitz::Array<cplx, 2> &A, const blitz::Array<cplx, 2> &B, blitz::Array<cplx, 2> &C);
void MatrixMatrixMultiply(const blitz::Array<double, 2> &A, const blitz::Array<double, 2> &B, blitz::Array<double, 2> &C);

//Performs the matrix-vector product w = A v
void MatrixVectorMultiply(const blitz::Array<cplx, 2> &A, const blitz::Array<cplx, 1> &v, blitz::Array<cplx, 1> &w);
void MatrixVectorMultiply(const blitz::Array<double, 2> &A, const blitz::Array<double, 1> &v, blitz::Array<double, 1> &w);

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

#endif


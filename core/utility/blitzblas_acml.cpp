
/*
 * This file should not be compiled separatly, it is included in blitzblas.cpp if
 * it is to be used
 */


//Use acml interface
namespace acml
{
#include <acml.h>
}
#define BLAS_NAME(name) acml::name


//Performs the matrix-matrix product C = A B
void MatrixMatrixMultiply(const Array<cplx, 2> &A, const Array<cplx, 2> &B, Array<cplx, 2> &C)
{
	/* We have c-style ordering, which corresponds to having all arrays transposed in 
	 * fortran-style ordering 
	 * and must therefore swich ordering and use transpose
	 * C := A B => C' := (B' A')
	*/

	//In blas world, C := C'
	int M = C.extent(1);
	int N = C.extent(0);
	int K = B.extent(0);

	int lda = A.stride(0);
	int ldb = B.stride(0);
	int ldc = C.stride(0);

	acml::doublecomplex alpha, beta;
	alpha.real = 1;
	alpha.imag = 0;
	beta.real = 0;
	beta.imag = 0;

	BLAS_NAME(zgemm)(
		'N',				// Transpose A
		'N',				// Transpose B
		M,				// Rows in C'
		N,				// Cols in C'
		K,				// Cols in A
		&alpha, 			// Scaling of A*B
		(acml::doublecomplex*)B.data(),	// Pointer to A
		ldb,				// Size of first dim of A
		(acml::doublecomplex*)A.data(),	// Pointer to B
		lda,				// Size of first dim of B
		&beta,				// Scaling of C
		(acml::doublecomplex*)C.data(),	// Pointer to C
		ldc				// Size of first dim of C
	);
}
		
void MatrixMatrixMultiply(const Array<double, 2> &A, const Array<double, 2> &B, Array<double, 2> &C)
{
	/* We have c-style ordering, but fortran-style ordering 
	 * and must therefore swich ordering and use transpose
	 * C := A B => C' := (B' A')
	*/

	//In blas world, C := C'
	int M = C.extent(1);
	int N = C.extent(0);
	int K = B.extent(0);

	int lda = A.stride(0);
	int ldb = B.stride(0);
	int ldc = C.stride(0);


	double alpha = 1;
	double beta = 0;

	BLAS_NAME(dgemm)(
		'N',				// Transpose A
		'N',				// Transpose B
		M,				// Rows in C'
		N,				// Cols in C'
		K,				// Cols in B'
		alpha, 				// Scaling of B' A'
		(double*)B.data(),		// Pointer to B'
		ldb,				// Size of first dim of B'
		(double*)A.data(),		// Pointer to A'
		lda,				// Size of first dim of A'
		beta,				// Scaling of C'
		(double*)C.data(),		// Pointer to C'
		ldc				// Size of first dim of C
	);
}

//Performs the matrix-vector product w = A v
void MatrixVectorMultiply(const Array<cplx, 2> &A, const Array<cplx, 1> &v, Array<cplx, 1> &w)
{
	//In lapack world, A := A', however, if we use 'T' to transpose A, all goes back to normal
	int M = A.extent(0);
	int N = A.extent(1);
	int lda = A.stride(0);

	int vStride = v.stride(0); 
	int wStride = w.stride(0); 

	acml::doublecomplex alpha, beta;
	alpha.real = 1;
	alpha.imag = 0;
	beta.real = 0;
	beta.imag = 0;

	BLAS_NAME(zgemv)(
		'T',					// Transpose A
		M,					// Rows of A 
		N, 					// Cols of A
		&alpha, 				// Scale of A*v
		(acml::doublecomplex*)A.data(), 				// Pointer to A
		lda, 					// Size of first dim of A
		(acml::doublecomplex*)v.data(),				// Pointer to input vector 
		vStride, 				// Stride of input vector
		&beta,	 				// Scale of w
		(acml::doublecomplex*)w.data(), 				// Pointer to output vector
		wStride					// Stride of output vector
	);
}


void MatrixVectorMultiply(const Array<double, 2> &A, const Array<double, 1> &v, Array<double, 1> &w)
{
	//In lapack world, A := A', however, if we use 'T' to transpose A, all goes back to normal
	int M = A.extent(0);
	int N = A.extent(1);
	int lda = A.stride(0);

	int vStride = v.stride(0);
	int wStride = w.stride(0);

	double alpha = 1;
	double beta = 0;

	BLAS_NAME(dgemv)(
		'T',		 			// Transpose A
		M,					// Rows of A 
		N, 					// Cols of A
		alpha, 					// Scale of A*v
		(double*)A.data(), 				// Pointer to A
		lda, 					// Size of first dim of A
		(double*)v.data(),				// Pointer to input vector 
		vStride, 				// Stride of input vector
		beta,	 				// Scale of w
		(double*)w.data(), 				// Pointer to output vector
		wStride					// Stride of output vector
	);		
}

//Performs elementwise multiplication of u and v
//There are no blas calls to do this, so we use the reference
//implementation
template<class T, int Rank>
void VectorElementMultiplyTemplate(const Array<T, Rank> &u, const Array<T, Rank> &v, Array<T, Rank> &w)
{
	w = u * v;
}

template<int Rank>
void VectorElementMultiply(const Array<double, Rank> &u, const Array<double, Rank> &v, Array<double, Rank> &w)
{
	VectorElementMultiplyTemplate(u, v, w);
}

template<int Rank>
void VectorElementMultiply(const Array<cplx, Rank> &u, const Array<cplx, Rank> &v, Array<cplx, Rank> &w)
{
	VectorElementMultiplyTemplate(u, v, w);
}


//Inner product of u and v
template<int Rank>
double VectorInnerProduct(const blitz::Array<double, Rank> &u, const blitz::Array<double, Rank> &v)
{
	if (u.size() != v.size())
	{
		cout << "Vector u and v is of different size: " << u.size() << " != " << v.size() << endl;
		throw std::runtime_error("invalid vector sizes for inner product");
	}
	if (!u.isStorageContiguous())
	{
		throw std::runtime_error("Vector u is not contiguous");
	}
	if (!v.isStorageContiguous())
	{
		throw std::runtime_error("Vector v is not contiguous");
	}
	int N = u.size();
	int uStride = 1; //u.stride(0);
	int vStride = 1; //v.stride(0);
	return BLAS_NAME(ddot)(N, (double*)u.data(), uStride, (double*)v.data(), vStride);
}

template<int Rank>
cplx VectorInnerProduct(const blitz::Array<cplx, Rank> &u, const blitz::Array<cplx, Rank> &v)
{
	if (u.size() != v.size())
	{
		cout << "Vector u and v is of different size: " << u.size() << " != " << v.size() << endl;
		throw std::runtime_error("invalid vector sizes for inner product");
	}
	if (!u.isStorageContiguous())
	{
		throw std::runtime_error("Vector u is not contiguous");
	}
	if (!v.isStorageContiguous())
	{
		throw std::runtime_error("Vector v is not contiguous");
	}
	int N = u.size();
	int uStride = 1; //u.stride(0);
	int vStride = 1; //v.stride(0);

	cplx result;
	acml::doublecomplex ret = BLAS_NAME(zdotc)(N, (acml::doublecomplex*)u.data(), uStride, (acml::doublecomplex*)v.data(), vStride);
	result = cplx(ret.real, ret.imag);
	return result;
}


//Explicitly instantiate template functions
template void VectorElementMultiply(const Array<cplx, 1> &u, const Array<cplx, 1> &v, Array<cplx, 1> &w);
template void VectorElementMultiply(const Array<cplx, 2> &u, const Array<cplx, 2> &v, Array<cplx, 2> &w);
template void VectorElementMultiply(const Array<cplx, 3> &u, const Array<cplx, 3> &v, Array<cplx, 3> &w);
template void VectorElementMultiply(const Array<cplx, 4> &u, const Array<cplx, 4> &v, Array<cplx, 4> &w);

template void VectorElementMultiply(const Array<double, 1> &u, const Array<double, 1> &v, Array<double, 1> &w);
template void VectorElementMultiply(const Array<double, 2> &u, const Array<double, 2> &v, Array<double, 2> &w);
template void VectorElementMultiply(const Array<double, 3> &u, const Array<double, 3> &v, Array<double, 3> &w);
template void VectorElementMultiply(const Array<double, 4> &u, const Array<double, 4> &v, Array<double, 4> &w);

template cplx VectorInnerProduct(const Array<cplx, 1> &u, const Array<cplx, 1> &v);
template cplx VectorInnerProduct(const Array<cplx, 2> &u, const Array<cplx, 2> &v);
template cplx VectorInnerProduct(const Array<cplx, 3> &u, const Array<cplx, 3> &v);
template cplx VectorInnerProduct(const Array<cplx, 4> &u, const Array<cplx, 4> &v);

template double VectorInnerProduct(const Array<double, 1> &u, const Array<double, 1> &v);
template double VectorInnerProduct(const Array<double, 2> &u, const Array<double, 2> &v);
template double VectorInnerProduct(const Array<double, 3> &u, const Array<double, 3> &v);
template double VectorInnerProduct(const Array<double, 4> &u, const Array<double, 4> &v);



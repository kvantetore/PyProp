
/*
 * This file should not be compiled separatly, it is included in blitzblas.cpp if
 * it is to be used
 */

//Performs the matrix-matrix product C = A B
template<class T>
void MatrixMatrixMultiplyTemplate(const Array<T, 2> &A, const Array<T, 2> &B, Array<T, 2> &C)
{
	using blitz::tensor::i;
	using blitz::tensor::j;
	using blitz::tensor::k;

	C = sum(A(i, k) * B(k, j), k);
}

void MatrixMatrixMultiply(const Array<cplx, 2> &A, const Array<cplx, 2> &B, Array<cplx, 2> &C)
{
	MatrixMatrixMultiplyTemplate(A, B, C);
}
		
void MatrixMatrixMultiply(const Array<double, 2> &A, const Array<double, 2> &B, Array<double, 2> &C)
{
	MatrixMatrixMultiplyTemplate(A, B, C);
}

//Performs the matrix-vector product w = A v
//
template<class T>
void MatrixVectorMultiplyTemplate(const Array<T, 2> &A, const Array<T, 1> &v, Array<T, 1> &w)
{
	using blitz::tensor::i;
	using blitz::tensor::j;
	w = sum(A(i, j) * v(j), j);
}

void MatrixVectorMultiply(const Array<cplx, 2> &A, const Array<cplx, 1> &v, Array<cplx, 1> &w)
{
	MatrixVectorMultiplyTemplate(A, v, w);
}

void MatrixVectorMultiply(const Array<double, 2> &A, const Array<double, 1> &v, Array<double, 1> &w)
{
	MatrixVectorMultiplyTemplate(A, v, w);
}

//Performs elementwise multiplication of u and v
template<class T, int Rank>
void VectorElementMultiplyTemplate(const Array<T, Rank> &u, const Array<T, Rank> &v, Array<T, Rank> &w)
{
	w = u * v;
	//cout << "hei\n";
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
	return sum(u * v);
}

template<int Rank>
cplx VectorInnerProduct(const blitz::Array<cplx, Rank> &u, const blitz::Array<cplx, Rank> &v)
{
	return sum(conj(u) * v);
}

//explicitly instantiate template functions
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



// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <utility/blitzblas.h>

// Using =======================================================================
using namespace boost::python;

void MatrixMatrixMultiply_cplx(blitz::Array<cplx, 2> A, blitz::Array<cplx, 2> B, blitz::Array<cplx, 2> C)
{
	MatrixMatrixMultiply(A, B, C);
}

void MatrixMatrixMultiply_double(blitz::Array<double, 2> A, blitz::Array<double, 2> B, blitz::Array<double, 2> C)
{
	MatrixMatrixMultiply(A, B, C);
}

void MatrixVectorMultiply_cplx(blitz::Array<cplx, 2> A, blitz::Array<cplx, 1> v, blitz::Array<cplx, 1> w)
{
	MatrixVectorMultiply(A, v, w);
}

void MatrixVectorMultiply_double(blitz::Array<double, 2> A, blitz::Array<double, 1> v, blitz::Array<double, 1> w)
{
	MatrixVectorMultiply(A, v, w);
}

void MatrixVectorMultiplyHermitianBanded_cplx(blitz::Array<cplx, 2> A, blitz::Array<cplx, 1> x, blitz::Array<cplx, 1> y,
	cplx alpha, cplx beta)
{
	MatrixVectorMultiplyHermitianBanded(A, x, y, alpha, beta);
}

/*
template<int Rank>
void VectorElementMultiply(blitz::Array<cplx, Rank> u, blitz::Array<cplx, Rank> v, blitz::Array<cplx, Rank> w);
template<int Rank>
void VectorElementMultiply(blitz::Array<double, Rank> u, blitz::Array<double, Rank> v, blitz::Array<double, Rank> w);

template<int Rank>
double VectorInnerProduct(blitz::Array<double, Rank> u, blitz::Array<double, Rank> v);
template<int Rank>
cplx VectorInnerProduct(blitz::Array<cplx, Rank> u, blitz::Array<cplx, Rank> v);
*/


// Module ======================================================================
void ExportBlitzBlas()
{
	def("MatrixMatrixMultiply", MatrixMatrixMultiply_double);
	def("MatrixVectorMultiply", MatrixVectorMultiply_double);
	def("MatrixVectorMultiplyHermitianBanded", MatrixVectorMultiplyHermitianBanded_cplx);
	//def("VectorElementMultiply_1", VectorElementMultiply<1>);
	//def("VectorInnerProduct_1", VectorInnerProduct<1>);
}


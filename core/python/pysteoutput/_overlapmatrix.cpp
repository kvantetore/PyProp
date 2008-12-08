
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/overlapmatrix.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

void (OverlapMatrix::*OverlapMatrix__MultiplyOverlapVectorconstblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1__)(const blitz::Array<std::complex<double>,1>&, blitz::Array<std::complex<double>,1>&)  = &OverlapMatrix::MultiplyOverlapVector;

void (OverlapMatrix::*OverlapMatrix__MultiplyOverlapVectorblitz__Array_std__complex_double__1__)(blitz::Array<std::complex<double>,1>&)  = &OverlapMatrix::MultiplyOverlapVector;

void (OverlapMatrix::*OverlapMatrix__MultiplyOverlapTensorconstblitz__Array_std__complex_double__3___blitz__Array_std__complex_double__3__)(const blitz::Array<std::complex<double>,3>&, blitz::Array<std::complex<double>,3>&)  = &OverlapMatrix::MultiplyOverlapTensor;

void (OverlapMatrix::*OverlapMatrix__MultiplyOverlapTensorblitz__Array_std__complex_double__3__)(blitz::Array<std::complex<double>,3>&)  = &OverlapMatrix::MultiplyOverlapTensor;


}// namespace 


// Module ======================================================================
void Export_python_overlapmatrix()
{
    scope* OverlapMatrix_scope = new scope(
    class_< OverlapMatrix, boost::noncopyable >("OverlapMatrix", init< int, int, const OverlapMatrixEvaluator& >())
        .def("GetOverlapHermitianUpper", &OverlapMatrix::GetOverlapHermitianUpper)
        .def("GetOverlapHermitianLower", &OverlapMatrix::GetOverlapHermitianLower)
        .def("GetOverlapCholeskyUpper", &OverlapMatrix::GetOverlapCholeskyUpper)
        .def("GetOverlapFullRow", &OverlapMatrix::GetOverlapFullRow)
        .def("GetOverlapFullCol", &OverlapMatrix::GetOverlapFullCol)
        .def("GetSuperDiagonals", &OverlapMatrix::GetSuperDiagonals)
        .def("GetBasisSize", &OverlapMatrix::GetBasisSize)
        .def("Setup", &OverlapMatrix::Setup)
        .def("MultiplyOverlapVector", OverlapMatrix__MultiplyOverlapVectorconstblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1__)
        .def("MultiplyOverlapVector", OverlapMatrix__MultiplyOverlapVectorblitz__Array_std__complex_double__1__)
        .def("SolveOverlapVector", &OverlapMatrix::SolveOverlapVector)
        .def("MultiplyOverlapTensor", OverlapMatrix__MultiplyOverlapTensorconstblitz__Array_std__complex_double__3___blitz__Array_std__complex_double__3__)
        .def("MultiplyOverlapTensor", OverlapMatrix__MultiplyOverlapTensorblitz__Array_std__complex_double__3__)
        .def("SolveOverlapTensor", &OverlapMatrix::SolveOverlapTensor)
        .def("MultiplySqrtOverlapVector", &OverlapMatrix::MultiplySqrtOverlapVector)
        .def("SolveSqrtOverlapVector", &OverlapMatrix::SolveSqrtOverlapVector)
        .def("MultiplySqrtOverlapTensor", &OverlapMatrix::MultiplySqrtOverlapTensor)
        .def("SolveSqrtOverlapTensor", &OverlapMatrix::SolveSqrtOverlapTensor)
    );
    register_ptr_to_python< boost::shared_ptr< OverlapMatrix > >();
    delete OverlapMatrix_scope;

}


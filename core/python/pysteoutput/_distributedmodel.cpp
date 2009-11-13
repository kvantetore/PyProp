
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <mpi/blitztranspose.h>
#include <mpi/distributedmodel.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

void (ArrayTranspose<1>::*ArrayTranspose_1___Transposeconstblitz__TinyVector_int_1__blitz__Array_std__complex_double__1__constblitz__Array_int_1__blitz__Array_std__complex_double__1__constblitz__Array_int_1_)(const blitz::TinyVector<int,1>, blitz::Array<std::complex<double>,1>, const blitz::Array<int,1>, blitz::Array<std::complex<double>,1>, const blitz::Array<int,1>)  = &ArrayTranspose<1>::Transpose;

void (ArrayTranspose<1>::*ArrayTranspose_1___Transposeconstblitz__TinyVector_int_1__constblitz__Array_int_1__blitz__Array_std__complex_double__1__int_blitz__Array_std__complex_double__1__int_int)(const blitz::TinyVector<int,1>, const blitz::Array<int,1>, blitz::Array<std::complex<double>,1>, int, blitz::Array<std::complex<double>,1>, int, int)  = &ArrayTranspose<1>::Transpose;

blitz::TinyVector<int,1> (ArrayTranspose<1>::*ArrayTranspose_1___CreateDistributedShapeconstblitz__TinyVector_int_1___constblitz__Array_int_1__)(const blitz::TinyVector<int,1>&, const blitz::Array<int,1>&)  = &ArrayTranspose<1>::CreateDistributedShape;

blitz::TinyVector<int,1> (ArrayTranspose<1>::*ArrayTranspose_1___CreateDistributedShapeconstblitz__TinyVector_int_1___constblitz__Array_int_1___constblitz__Array_int_1__)(const blitz::TinyVector<int,1>&, const blitz::Array<int,1>&, const blitz::Array<int,1>&)  = &ArrayTranspose<1>::CreateDistributedShape;

blitz::TinyVector<int,1> (ArrayTranspose<1>::*ArrayTranspose_1___CreatePaddedShapeconstblitz__TinyVector_int_1___constblitz__Array_int_1__)(const blitz::TinyVector<int,1>&, const blitz::Array<int,1>&)  = &ArrayTranspose<1>::CreatePaddedShape;

int (ArrayTranspose<1>::*ArrayTranspose_1___CreatePaddedShapeint_int)(int, int)  = &ArrayTranspose<1>::CreatePaddedShape;

int (ArrayTranspose<1>::*ArrayTranspose_1___CreateDistributedShapeint_int)(int, int)  = &ArrayTranspose<1>::CreateDistributedShape;

int (ArrayTranspose<1>::*ArrayTranspose_1___CreateDistributedShapeint_int_int)(int, int, int)  = &ArrayTranspose<1>::CreateDistributedShape;

int (ArrayTranspose<1>::*ArrayTranspose_1___GetLocalStartIndexint_int_int)(int, int, int)  = &ArrayTranspose<1>::GetLocalStartIndex;

blitz::Range (ArrayTranspose<1>::*ArrayTranspose_1___GetLocalRangeint_int_int)(int, int, int)  = &ArrayTranspose<1>::GetLocalRange;

int (ArrayTranspose<1>::*ArrayTranspose_1___GetLocalStartIndexint_int)(int, int)  = &ArrayTranspose<1>::GetLocalStartIndex;

blitz::Range (ArrayTranspose<1>::*ArrayTranspose_1___GetLocalRangeint_int)(int, int)  = &ArrayTranspose<1>::GetLocalRange;

double (ArrayTranspose<1>::*ArrayTranspose_1___GetGlobalSumdouble_int)(double, int)  = &ArrayTranspose<1>::GetGlobalSum;

std::complex<double> (ArrayTranspose<1>::*ArrayTranspose_1___GetGlobalSumstd__complex_double__int)(std::complex<double>, int)  = &ArrayTranspose<1>::GetGlobalSum;

int (ArrayTranspose<1>::*ArrayTranspose_1___GetGlobalSumint_int)(int, int)  = &ArrayTranspose<1>::GetGlobalSum;

void (ArrayTranspose<2>::*ArrayTranspose_2___Transposeconstblitz__TinyVector_int_2__blitz__Array_std__complex_double__2__constblitz__Array_int_1__blitz__Array_std__complex_double__2__constblitz__Array_int_1_)(const blitz::TinyVector<int,2>, blitz::Array<std::complex<double>,2>, const blitz::Array<int,1>, blitz::Array<std::complex<double>,2>, const blitz::Array<int,1>)  = &ArrayTranspose<2>::Transpose;

void (ArrayTranspose<2>::*ArrayTranspose_2___Transposeconstblitz__TinyVector_int_2__constblitz__Array_int_1__blitz__Array_std__complex_double__2__int_blitz__Array_std__complex_double__2__int_int)(const blitz::TinyVector<int,2>, const blitz::Array<int,1>, blitz::Array<std::complex<double>,2>, int, blitz::Array<std::complex<double>,2>, int, int)  = &ArrayTranspose<2>::Transpose;

blitz::TinyVector<int,2> (ArrayTranspose<2>::*ArrayTranspose_2___CreateDistributedShapeconstblitz__TinyVector_int_2___constblitz__Array_int_1__)(const blitz::TinyVector<int,2>&, const blitz::Array<int,1>&)  = &ArrayTranspose<2>::CreateDistributedShape;

blitz::TinyVector<int,2> (ArrayTranspose<2>::*ArrayTranspose_2___CreateDistributedShapeconstblitz__TinyVector_int_2___constblitz__Array_int_1___constblitz__Array_int_1__)(const blitz::TinyVector<int,2>&, const blitz::Array<int,1>&, const blitz::Array<int,1>&)  = &ArrayTranspose<2>::CreateDistributedShape;

blitz::TinyVector<int,2> (ArrayTranspose<2>::*ArrayTranspose_2___CreatePaddedShapeconstblitz__TinyVector_int_2___constblitz__Array_int_1__)(const blitz::TinyVector<int,2>&, const blitz::Array<int,1>&)  = &ArrayTranspose<2>::CreatePaddedShape;

int (ArrayTranspose<2>::*ArrayTranspose_2___CreatePaddedShapeint_int)(int, int)  = &ArrayTranspose<2>::CreatePaddedShape;

int (ArrayTranspose<2>::*ArrayTranspose_2___CreateDistributedShapeint_int)(int, int)  = &ArrayTranspose<2>::CreateDistributedShape;

int (ArrayTranspose<2>::*ArrayTranspose_2___CreateDistributedShapeint_int_int)(int, int, int)  = &ArrayTranspose<2>::CreateDistributedShape;

int (ArrayTranspose<2>::*ArrayTranspose_2___GetLocalStartIndexint_int_int)(int, int, int)  = &ArrayTranspose<2>::GetLocalStartIndex;

blitz::Range (ArrayTranspose<2>::*ArrayTranspose_2___GetLocalRangeint_int_int)(int, int, int)  = &ArrayTranspose<2>::GetLocalRange;

int (ArrayTranspose<2>::*ArrayTranspose_2___GetLocalStartIndexint_int)(int, int)  = &ArrayTranspose<2>::GetLocalStartIndex;

blitz::Range (ArrayTranspose<2>::*ArrayTranspose_2___GetLocalRangeint_int)(int, int)  = &ArrayTranspose<2>::GetLocalRange;

double (ArrayTranspose<2>::*ArrayTranspose_2___GetGlobalSumdouble_int)(double, int)  = &ArrayTranspose<2>::GetGlobalSum;

std::complex<double> (ArrayTranspose<2>::*ArrayTranspose_2___GetGlobalSumstd__complex_double__int)(std::complex<double>, int)  = &ArrayTranspose<2>::GetGlobalSum;

int (ArrayTranspose<2>::*ArrayTranspose_2___GetGlobalSumint_int)(int, int)  = &ArrayTranspose<2>::GetGlobalSum;

void (ArrayTranspose<3>::*ArrayTranspose_3___Transposeconstblitz__TinyVector_int_3__blitz__Array_std__complex_double__3__constblitz__Array_int_1__blitz__Array_std__complex_double__3__constblitz__Array_int_1_)(const blitz::TinyVector<int,3>, blitz::Array<std::complex<double>,3>, const blitz::Array<int,1>, blitz::Array<std::complex<double>,3>, const blitz::Array<int,1>)  = &ArrayTranspose<3>::Transpose;

void (ArrayTranspose<3>::*ArrayTranspose_3___Transposeconstblitz__TinyVector_int_3__constblitz__Array_int_1__blitz__Array_std__complex_double__3__int_blitz__Array_std__complex_double__3__int_int)(const blitz::TinyVector<int,3>, const blitz::Array<int,1>, blitz::Array<std::complex<double>,3>, int, blitz::Array<std::complex<double>,3>, int, int)  = &ArrayTranspose<3>::Transpose;

blitz::TinyVector<int,3> (ArrayTranspose<3>::*ArrayTranspose_3___CreateDistributedShapeconstblitz__TinyVector_int_3___constblitz__Array_int_1__)(const blitz::TinyVector<int,3>&, const blitz::Array<int,1>&)  = &ArrayTranspose<3>::CreateDistributedShape;

blitz::TinyVector<int,3> (ArrayTranspose<3>::*ArrayTranspose_3___CreateDistributedShapeconstblitz__TinyVector_int_3___constblitz__Array_int_1___constblitz__Array_int_1__)(const blitz::TinyVector<int,3>&, const blitz::Array<int,1>&, const blitz::Array<int,1>&)  = &ArrayTranspose<3>::CreateDistributedShape;

blitz::TinyVector<int,3> (ArrayTranspose<3>::*ArrayTranspose_3___CreatePaddedShapeconstblitz__TinyVector_int_3___constblitz__Array_int_1__)(const blitz::TinyVector<int,3>&, const blitz::Array<int,1>&)  = &ArrayTranspose<3>::CreatePaddedShape;

int (ArrayTranspose<3>::*ArrayTranspose_3___CreatePaddedShapeint_int)(int, int)  = &ArrayTranspose<3>::CreatePaddedShape;

int (ArrayTranspose<3>::*ArrayTranspose_3___CreateDistributedShapeint_int)(int, int)  = &ArrayTranspose<3>::CreateDistributedShape;

int (ArrayTranspose<3>::*ArrayTranspose_3___CreateDistributedShapeint_int_int)(int, int, int)  = &ArrayTranspose<3>::CreateDistributedShape;

int (ArrayTranspose<3>::*ArrayTranspose_3___GetLocalStartIndexint_int_int)(int, int, int)  = &ArrayTranspose<3>::GetLocalStartIndex;

blitz::Range (ArrayTranspose<3>::*ArrayTranspose_3___GetLocalRangeint_int_int)(int, int, int)  = &ArrayTranspose<3>::GetLocalRange;

int (ArrayTranspose<3>::*ArrayTranspose_3___GetLocalStartIndexint_int)(int, int)  = &ArrayTranspose<3>::GetLocalStartIndex;

blitz::Range (ArrayTranspose<3>::*ArrayTranspose_3___GetLocalRangeint_int)(int, int)  = &ArrayTranspose<3>::GetLocalRange;

double (ArrayTranspose<3>::*ArrayTranspose_3___GetGlobalSumdouble_int)(double, int)  = &ArrayTranspose<3>::GetGlobalSum;

std::complex<double> (ArrayTranspose<3>::*ArrayTranspose_3___GetGlobalSumstd__complex_double__int)(std::complex<double>, int)  = &ArrayTranspose<3>::GetGlobalSum;

int (ArrayTranspose<3>::*ArrayTranspose_3___GetGlobalSumint_int)(int, int)  = &ArrayTranspose<3>::GetGlobalSum;

void (ArrayTranspose<4>::*ArrayTranspose_4___Transposeconstblitz__TinyVector_int_4__blitz__Array_std__complex_double__4__constblitz__Array_int_1__blitz__Array_std__complex_double__4__constblitz__Array_int_1_)(const blitz::TinyVector<int,4>, blitz::Array<std::complex<double>,4>, const blitz::Array<int,1>, blitz::Array<std::complex<double>,4>, const blitz::Array<int,1>)  = &ArrayTranspose<4>::Transpose;

void (ArrayTranspose<4>::*ArrayTranspose_4___Transposeconstblitz__TinyVector_int_4__constblitz__Array_int_1__blitz__Array_std__complex_double__4__int_blitz__Array_std__complex_double__4__int_int)(const blitz::TinyVector<int,4>, const blitz::Array<int,1>, blitz::Array<std::complex<double>,4>, int, blitz::Array<std::complex<double>,4>, int, int)  = &ArrayTranspose<4>::Transpose;

blitz::TinyVector<int,4> (ArrayTranspose<4>::*ArrayTranspose_4___CreateDistributedShapeconstblitz__TinyVector_int_4___constblitz__Array_int_1__)(const blitz::TinyVector<int,4>&, const blitz::Array<int,1>&)  = &ArrayTranspose<4>::CreateDistributedShape;

blitz::TinyVector<int,4> (ArrayTranspose<4>::*ArrayTranspose_4___CreateDistributedShapeconstblitz__TinyVector_int_4___constblitz__Array_int_1___constblitz__Array_int_1__)(const blitz::TinyVector<int,4>&, const blitz::Array<int,1>&, const blitz::Array<int,1>&)  = &ArrayTranspose<4>::CreateDistributedShape;

blitz::TinyVector<int,4> (ArrayTranspose<4>::*ArrayTranspose_4___CreatePaddedShapeconstblitz__TinyVector_int_4___constblitz__Array_int_1__)(const blitz::TinyVector<int,4>&, const blitz::Array<int,1>&)  = &ArrayTranspose<4>::CreatePaddedShape;

int (ArrayTranspose<4>::*ArrayTranspose_4___CreatePaddedShapeint_int)(int, int)  = &ArrayTranspose<4>::CreatePaddedShape;

int (ArrayTranspose<4>::*ArrayTranspose_4___CreateDistributedShapeint_int)(int, int)  = &ArrayTranspose<4>::CreateDistributedShape;

int (ArrayTranspose<4>::*ArrayTranspose_4___CreateDistributedShapeint_int_int)(int, int, int)  = &ArrayTranspose<4>::CreateDistributedShape;

int (ArrayTranspose<4>::*ArrayTranspose_4___GetLocalStartIndexint_int_int)(int, int, int)  = &ArrayTranspose<4>::GetLocalStartIndex;

blitz::Range (ArrayTranspose<4>::*ArrayTranspose_4___GetLocalRangeint_int_int)(int, int, int)  = &ArrayTranspose<4>::GetLocalRange;

int (ArrayTranspose<4>::*ArrayTranspose_4___GetLocalStartIndexint_int)(int, int)  = &ArrayTranspose<4>::GetLocalStartIndex;

blitz::Range (ArrayTranspose<4>::*ArrayTranspose_4___GetLocalRangeint_int)(int, int)  = &ArrayTranspose<4>::GetLocalRange;

double (ArrayTranspose<4>::*ArrayTranspose_4___GetGlobalSumdouble_int)(double, int)  = &ArrayTranspose<4>::GetGlobalSum;

std::complex<double> (ArrayTranspose<4>::*ArrayTranspose_4___GetGlobalSumstd__complex_double__int)(std::complex<double>, int)  = &ArrayTranspose<4>::GetGlobalSum;

int (ArrayTranspose<4>::*ArrayTranspose_4___GetGlobalSumint_int)(int, int)  = &ArrayTranspose<4>::GetGlobalSum;

int (DistributedModel<1>::*DistributedModel_1___GetLocalStartIndexint_int)(int, int)  = &DistributedModel<1>::GetLocalStartIndex;

blitz::Range (DistributedModel<1>::*DistributedModel_1___GetLocalIndexRangeint_int)(int, int)  = &DistributedModel<1>::GetLocalIndexRange;

int (DistributedModel<1>::*DistributedModel_1___GetLocalStartIndexint_int_int)(int, int, int)  = &DistributedModel<1>::GetLocalStartIndex;

blitz::Range (DistributedModel<1>::*DistributedModel_1___GetLocalIndexRangeint_int_int)(int, int, int)  = &DistributedModel<1>::GetLocalIndexRange;

double (DistributedModel<1>::*DistributedModel_1___GetGlobalSumdouble)(double)  = &DistributedModel<1>::GetGlobalSum;

std::complex<double> (DistributedModel<1>::*DistributedModel_1___GetGlobalSumstd__complex_double_)(std::complex<double>)  = &DistributedModel<1>::GetGlobalSum;

void (DistributedModel<1>::*DistributedModel_1___GetGlobalSumblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1__)(blitz::Array<std::complex<double>,1>&, blitz::Array<std::complex<double>,1>&)  = &DistributedModel<1>::GetGlobalSum;

int (DistributedModel<2>::*DistributedModel_2___GetLocalStartIndexint_int)(int, int)  = &DistributedModel<2>::GetLocalStartIndex;

blitz::Range (DistributedModel<2>::*DistributedModel_2___GetLocalIndexRangeint_int)(int, int)  = &DistributedModel<2>::GetLocalIndexRange;

int (DistributedModel<2>::*DistributedModel_2___GetLocalStartIndexint_int_int)(int, int, int)  = &DistributedModel<2>::GetLocalStartIndex;

blitz::Range (DistributedModel<2>::*DistributedModel_2___GetLocalIndexRangeint_int_int)(int, int, int)  = &DistributedModel<2>::GetLocalIndexRange;

double (DistributedModel<2>::*DistributedModel_2___GetGlobalSumdouble)(double)  = &DistributedModel<2>::GetGlobalSum;

std::complex<double> (DistributedModel<2>::*DistributedModel_2___GetGlobalSumstd__complex_double_)(std::complex<double>)  = &DistributedModel<2>::GetGlobalSum;

void (DistributedModel<2>::*DistributedModel_2___GetGlobalSumblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1__)(blitz::Array<std::complex<double>,1>&, blitz::Array<std::complex<double>,1>&)  = &DistributedModel<2>::GetGlobalSum;

int (DistributedModel<3>::*DistributedModel_3___GetLocalStartIndexint_int)(int, int)  = &DistributedModel<3>::GetLocalStartIndex;

blitz::Range (DistributedModel<3>::*DistributedModel_3___GetLocalIndexRangeint_int)(int, int)  = &DistributedModel<3>::GetLocalIndexRange;

int (DistributedModel<3>::*DistributedModel_3___GetLocalStartIndexint_int_int)(int, int, int)  = &DistributedModel<3>::GetLocalStartIndex;

blitz::Range (DistributedModel<3>::*DistributedModel_3___GetLocalIndexRangeint_int_int)(int, int, int)  = &DistributedModel<3>::GetLocalIndexRange;

double (DistributedModel<3>::*DistributedModel_3___GetGlobalSumdouble)(double)  = &DistributedModel<3>::GetGlobalSum;

std::complex<double> (DistributedModel<3>::*DistributedModel_3___GetGlobalSumstd__complex_double_)(std::complex<double>)  = &DistributedModel<3>::GetGlobalSum;

void (DistributedModel<3>::*DistributedModel_3___GetGlobalSumblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1__)(blitz::Array<std::complex<double>,1>&, blitz::Array<std::complex<double>,1>&)  = &DistributedModel<3>::GetGlobalSum;

int (DistributedModel<4>::*DistributedModel_4___GetLocalStartIndexint_int)(int, int)  = &DistributedModel<4>::GetLocalStartIndex;

blitz::Range (DistributedModel<4>::*DistributedModel_4___GetLocalIndexRangeint_int)(int, int)  = &DistributedModel<4>::GetLocalIndexRange;

int (DistributedModel<4>::*DistributedModel_4___GetLocalStartIndexint_int_int)(int, int, int)  = &DistributedModel<4>::GetLocalStartIndex;

blitz::Range (DistributedModel<4>::*DistributedModel_4___GetLocalIndexRangeint_int_int)(int, int, int)  = &DistributedModel<4>::GetLocalIndexRange;

double (DistributedModel<4>::*DistributedModel_4___GetGlobalSumdouble)(double)  = &DistributedModel<4>::GetGlobalSum;

std::complex<double> (DistributedModel<4>::*DistributedModel_4___GetGlobalSumstd__complex_double_)(std::complex<double>)  = &DistributedModel<4>::GetGlobalSum;

void (DistributedModel<4>::*DistributedModel_4___GetGlobalSumblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1__)(blitz::Array<std::complex<double>,1>&, blitz::Array<std::complex<double>,1>&)  = &DistributedModel<4>::GetGlobalSum;


}// namespace 


// Module ======================================================================
void Export_python_distributedmodel()
{
    scope* ArrayTranspose_1_scope = new scope(
    class_< ArrayTranspose<1>, boost::noncopyable >("ArrayTranspose_1", init< int >())
        .def("GetProcGridShape", &ArrayTranspose<1>::GetProcGridShape, return_value_policy< copy_const_reference >())
        .def("Transpose", ArrayTranspose_1___Transposeconstblitz__TinyVector_int_1__blitz__Array_std__complex_double__1__constblitz__Array_int_1__blitz__Array_std__complex_double__1__constblitz__Array_int_1_)
        .def("Transpose", ArrayTranspose_1___Transposeconstblitz__TinyVector_int_1__constblitz__Array_int_1__blitz__Array_std__complex_double__1__int_blitz__Array_std__complex_double__1__int_int)
        .def("CreateDistributedShape", ArrayTranspose_1___CreateDistributedShapeconstblitz__TinyVector_int_1___constblitz__Array_int_1__)
        .def("CreateDistributedShape", ArrayTranspose_1___CreateDistributedShapeconstblitz__TinyVector_int_1___constblitz__Array_int_1___constblitz__Array_int_1__)
        .def("CreatePaddedShape", ArrayTranspose_1___CreatePaddedShapeconstblitz__TinyVector_int_1___constblitz__Array_int_1__)
        .def("CreatePaddedShape", ArrayTranspose_1___CreatePaddedShapeint_int)
        .def("CreateDistributedShape", ArrayTranspose_1___CreateDistributedShapeint_int)
        .def("CreateDistributedShape", ArrayTranspose_1___CreateDistributedShapeint_int_int)
        .def("GetLocalStartIndex", ArrayTranspose_1___GetLocalStartIndexint_int_int)
        .def("GetLocalRange", ArrayTranspose_1___GetLocalRangeint_int_int)
        .def("GetLocalStartIndex", ArrayTranspose_1___GetLocalStartIndexint_int)
        .def("GetLocalRange", ArrayTranspose_1___GetLocalRangeint_int)
        .def("CreateReferenceData", &ArrayTranspose<1>::CreateReferenceData)
        .def("GetGlobalSum", ArrayTranspose_1___GetGlobalSumdouble_int)
        .def("GetGlobalSum", ArrayTranspose_1___GetGlobalSumstd__complex_double__int)
        .def("GetGlobalSum", ArrayTranspose_1___GetGlobalSumint_int)
        .def("GetGroupComm", &ArrayTranspose<1>::GetGroupComm)
    );
    register_ptr_to_python< boost::shared_ptr< ArrayTranspose<1> > >();
    delete ArrayTranspose_1_scope;

    scope* ArrayTranspose_2_scope = new scope(
    class_< ArrayTranspose<2>, boost::noncopyable >("ArrayTranspose_2", init< int >())
        .def("GetProcGridShape", &ArrayTranspose<2>::GetProcGridShape, return_value_policy< copy_const_reference >())
        .def("Transpose", ArrayTranspose_2___Transposeconstblitz__TinyVector_int_2__blitz__Array_std__complex_double__2__constblitz__Array_int_1__blitz__Array_std__complex_double__2__constblitz__Array_int_1_)
        .def("Transpose", ArrayTranspose_2___Transposeconstblitz__TinyVector_int_2__constblitz__Array_int_1__blitz__Array_std__complex_double__2__int_blitz__Array_std__complex_double__2__int_int)
        .def("CreateDistributedShape", ArrayTranspose_2___CreateDistributedShapeconstblitz__TinyVector_int_2___constblitz__Array_int_1__)
        .def("CreateDistributedShape", ArrayTranspose_2___CreateDistributedShapeconstblitz__TinyVector_int_2___constblitz__Array_int_1___constblitz__Array_int_1__)
        .def("CreatePaddedShape", ArrayTranspose_2___CreatePaddedShapeconstblitz__TinyVector_int_2___constblitz__Array_int_1__)
        .def("CreatePaddedShape", ArrayTranspose_2___CreatePaddedShapeint_int)
        .def("CreateDistributedShape", ArrayTranspose_2___CreateDistributedShapeint_int)
        .def("CreateDistributedShape", ArrayTranspose_2___CreateDistributedShapeint_int_int)
        .def("GetLocalStartIndex", ArrayTranspose_2___GetLocalStartIndexint_int_int)
        .def("GetLocalRange", ArrayTranspose_2___GetLocalRangeint_int_int)
        .def("GetLocalStartIndex", ArrayTranspose_2___GetLocalStartIndexint_int)
        .def("GetLocalRange", ArrayTranspose_2___GetLocalRangeint_int)
        .def("CreateReferenceData", &ArrayTranspose<2>::CreateReferenceData)
        .def("GetGlobalSum", ArrayTranspose_2___GetGlobalSumdouble_int)
        .def("GetGlobalSum", ArrayTranspose_2___GetGlobalSumstd__complex_double__int)
        .def("GetGlobalSum", ArrayTranspose_2___GetGlobalSumint_int)
        .def("GetGroupComm", &ArrayTranspose<2>::GetGroupComm)
    );
    register_ptr_to_python< boost::shared_ptr< ArrayTranspose<2> > >();
    delete ArrayTranspose_2_scope;

    scope* ArrayTranspose_3_scope = new scope(
    class_< ArrayTranspose<3>, boost::noncopyable >("ArrayTranspose_3", init< int >())
        .def("GetProcGridShape", &ArrayTranspose<3>::GetProcGridShape, return_value_policy< copy_const_reference >())
        .def("Transpose", ArrayTranspose_3___Transposeconstblitz__TinyVector_int_3__blitz__Array_std__complex_double__3__constblitz__Array_int_1__blitz__Array_std__complex_double__3__constblitz__Array_int_1_)
        .def("Transpose", ArrayTranspose_3___Transposeconstblitz__TinyVector_int_3__constblitz__Array_int_1__blitz__Array_std__complex_double__3__int_blitz__Array_std__complex_double__3__int_int)
        .def("CreateDistributedShape", ArrayTranspose_3___CreateDistributedShapeconstblitz__TinyVector_int_3___constblitz__Array_int_1__)
        .def("CreateDistributedShape", ArrayTranspose_3___CreateDistributedShapeconstblitz__TinyVector_int_3___constblitz__Array_int_1___constblitz__Array_int_1__)
        .def("CreatePaddedShape", ArrayTranspose_3___CreatePaddedShapeconstblitz__TinyVector_int_3___constblitz__Array_int_1__)
        .def("CreatePaddedShape", ArrayTranspose_3___CreatePaddedShapeint_int)
        .def("CreateDistributedShape", ArrayTranspose_3___CreateDistributedShapeint_int)
        .def("CreateDistributedShape", ArrayTranspose_3___CreateDistributedShapeint_int_int)
        .def("GetLocalStartIndex", ArrayTranspose_3___GetLocalStartIndexint_int_int)
        .def("GetLocalRange", ArrayTranspose_3___GetLocalRangeint_int_int)
        .def("GetLocalStartIndex", ArrayTranspose_3___GetLocalStartIndexint_int)
        .def("GetLocalRange", ArrayTranspose_3___GetLocalRangeint_int)
        .def("CreateReferenceData", &ArrayTranspose<3>::CreateReferenceData)
        .def("GetGlobalSum", ArrayTranspose_3___GetGlobalSumdouble_int)
        .def("GetGlobalSum", ArrayTranspose_3___GetGlobalSumstd__complex_double__int)
        .def("GetGlobalSum", ArrayTranspose_3___GetGlobalSumint_int)
        .def("GetGroupComm", &ArrayTranspose<3>::GetGroupComm)
    );
    register_ptr_to_python< boost::shared_ptr< ArrayTranspose<3> > >();
    delete ArrayTranspose_3_scope;

    scope* ArrayTranspose_4_scope = new scope(
    class_< ArrayTranspose<4>, boost::noncopyable >("ArrayTranspose_4", init< int >())
        .def("GetProcGridShape", &ArrayTranspose<4>::GetProcGridShape, return_value_policy< copy_const_reference >())
        .def("Transpose", ArrayTranspose_4___Transposeconstblitz__TinyVector_int_4__blitz__Array_std__complex_double__4__constblitz__Array_int_1__blitz__Array_std__complex_double__4__constblitz__Array_int_1_)
        .def("Transpose", ArrayTranspose_4___Transposeconstblitz__TinyVector_int_4__constblitz__Array_int_1__blitz__Array_std__complex_double__4__int_blitz__Array_std__complex_double__4__int_int)
        .def("CreateDistributedShape", ArrayTranspose_4___CreateDistributedShapeconstblitz__TinyVector_int_4___constblitz__Array_int_1__)
        .def("CreateDistributedShape", ArrayTranspose_4___CreateDistributedShapeconstblitz__TinyVector_int_4___constblitz__Array_int_1___constblitz__Array_int_1__)
        .def("CreatePaddedShape", ArrayTranspose_4___CreatePaddedShapeconstblitz__TinyVector_int_4___constblitz__Array_int_1__)
        .def("CreatePaddedShape", ArrayTranspose_4___CreatePaddedShapeint_int)
        .def("CreateDistributedShape", ArrayTranspose_4___CreateDistributedShapeint_int)
        .def("CreateDistributedShape", ArrayTranspose_4___CreateDistributedShapeint_int_int)
        .def("GetLocalStartIndex", ArrayTranspose_4___GetLocalStartIndexint_int_int)
        .def("GetLocalRange", ArrayTranspose_4___GetLocalRangeint_int_int)
        .def("GetLocalStartIndex", ArrayTranspose_4___GetLocalStartIndexint_int)
        .def("GetLocalRange", ArrayTranspose_4___GetLocalRangeint_int)
        .def("CreateReferenceData", &ArrayTranspose<4>::CreateReferenceData)
        .def("GetGlobalSum", ArrayTranspose_4___GetGlobalSumdouble_int)
        .def("GetGlobalSum", ArrayTranspose_4___GetGlobalSumstd__complex_double__int)
        .def("GetGlobalSum", ArrayTranspose_4___GetGlobalSumint_int)
        .def("GetGroupComm", &ArrayTranspose<4>::GetGroupComm)
    );
    register_ptr_to_python< boost::shared_ptr< ArrayTranspose<4> > >();
    delete ArrayTranspose_4_scope;

    scope* DistributedModel_1_scope = new scope(
    class_< DistributedModel<1> >("DistributedModel_1", init<  >())
        .def(init< boost::shared_ptr<Distribution> >())
        .def(init< const DistributedModel<1>& >())
        .def_readwrite("ProcId", &DistributedModel<1>::ProcId)
        .def_readwrite("ProcCount", &DistributedModel<1>::ProcCount)
        .def("CreateSubDistributedModel", &DistributedModel<1>::CreateSubDistributedModel)
        .def("IsSingleProc", &DistributedModel<1>::IsSingleProc)
        .def("IsFirstProc", &DistributedModel<1>::IsFirstProc)
        .def("IsLastProc", &DistributedModel<1>::IsLastProc)
        .def("GetDistribution", &DistributedModel<1>::GetDistribution)
        .def("SetDistribution", &DistributedModel<1>::SetDistribution)
        .def("GetTranspose", &DistributedModel<1>::GetTranspose)
        .def("ForceSingleProc", &DistributedModel<1>::ForceSingleProc)
        .def("SetupMPI", &DistributedModel<1>::SetupMPI)
        .def("CreateInitialShape", &DistributedModel<1>::CreateInitialShape)
        .def("GetGlobalShape", &DistributedModel<1>::GetGlobalShape)
        .def("GetLocalStartIndex", DistributedModel_1___GetLocalStartIndexint_int)
        .def("GetLocalIndexRange", DistributedModel_1___GetLocalIndexRangeint_int)
        .def("ApplyConfigSection", &DistributedModel<1>::ApplyConfigSection)
        .def("GetLocalStartIndex", DistributedModel_1___GetLocalStartIndexint_int_int)
        .def("GetLocalIndexRange", DistributedModel_1___GetLocalIndexRangeint_int_int)
        .def("GetGlobalSum", DistributedModel_1___GetGlobalSumdouble)
        .def("GetGlobalSum", DistributedModel_1___GetGlobalSumstd__complex_double_)
        .def("GetGlobalSum", DistributedModel_1___GetGlobalSumblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1__)
        .def("GlobalBarrier", &DistributedModel<1>::GlobalBarrier)
        .def("IsDistributedRank", &DistributedModel<1>::IsDistributedRank)
        .def("ChangeDistribution", &DistributedModel<1>::ChangeDistribution)
        .def("InitMPI", &DistributedModel<1>::InitMPI)
        .def("FinalizeMPI", &DistributedModel<1>::FinalizeMPI)
        .staticmethod("InitMPI")
        .staticmethod("FinalizeMPI")
        .staticmethod("ForceSingleProc")
    );
    register_ptr_to_python< boost::shared_ptr< DistributedModel<1> > >();
    delete DistributedModel_1_scope;

    scope* DistributedModel_2_scope = new scope(
    class_< DistributedModel<2> >("DistributedModel_2", init<  >())
        .def(init< boost::shared_ptr<Distribution> >())
        .def(init< const DistributedModel<2>& >())
        .def_readwrite("ProcId", &DistributedModel<2>::ProcId)
        .def_readwrite("ProcCount", &DistributedModel<2>::ProcCount)
        .def("CreateSubDistributedModel", &DistributedModel<2>::CreateSubDistributedModel)
        .def("IsSingleProc", &DistributedModel<2>::IsSingleProc)
        .def("IsFirstProc", &DistributedModel<2>::IsFirstProc)
        .def("IsLastProc", &DistributedModel<2>::IsLastProc)
        .def("GetDistribution", &DistributedModel<2>::GetDistribution)
        .def("SetDistribution", &DistributedModel<2>::SetDistribution)
        .def("GetTranspose", &DistributedModel<2>::GetTranspose)
        .def("ForceSingleProc", &DistributedModel<2>::ForceSingleProc)
        .def("SetupMPI", &DistributedModel<2>::SetupMPI)
        .def("CreateInitialShape", &DistributedModel<2>::CreateInitialShape)
        .def("GetGlobalShape", &DistributedModel<2>::GetGlobalShape)
        .def("GetLocalStartIndex", DistributedModel_2___GetLocalStartIndexint_int)
        .def("GetLocalIndexRange", DistributedModel_2___GetLocalIndexRangeint_int)
        .def("ApplyConfigSection", &DistributedModel<2>::ApplyConfigSection)
        .def("GetLocalStartIndex", DistributedModel_2___GetLocalStartIndexint_int_int)
        .def("GetLocalIndexRange", DistributedModel_2___GetLocalIndexRangeint_int_int)
        .def("GetGlobalSum", DistributedModel_2___GetGlobalSumdouble)
        .def("GetGlobalSum", DistributedModel_2___GetGlobalSumstd__complex_double_)
        .def("GetGlobalSum", DistributedModel_2___GetGlobalSumblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1__)
        .def("GlobalBarrier", &DistributedModel<2>::GlobalBarrier)
        .def("IsDistributedRank", &DistributedModel<2>::IsDistributedRank)
        .def("ChangeDistribution", &DistributedModel<2>::ChangeDistribution)
        .def("InitMPI", &DistributedModel<2>::InitMPI)
        .def("FinalizeMPI", &DistributedModel<2>::FinalizeMPI)
        .staticmethod("InitMPI")
        .staticmethod("FinalizeMPI")
        .staticmethod("ForceSingleProc")
    );
    register_ptr_to_python< boost::shared_ptr< DistributedModel<2> > >();
    delete DistributedModel_2_scope;

    scope* DistributedModel_3_scope = new scope(
    class_< DistributedModel<3> >("DistributedModel_3", init<  >())
        .def(init< boost::shared_ptr<Distribution> >())
        .def(init< const DistributedModel<3>& >())
        .def_readwrite("ProcId", &DistributedModel<3>::ProcId)
        .def_readwrite("ProcCount", &DistributedModel<3>::ProcCount)
        .def("CreateSubDistributedModel", &DistributedModel<3>::CreateSubDistributedModel)
        .def("IsSingleProc", &DistributedModel<3>::IsSingleProc)
        .def("IsFirstProc", &DistributedModel<3>::IsFirstProc)
        .def("IsLastProc", &DistributedModel<3>::IsLastProc)
        .def("GetDistribution", &DistributedModel<3>::GetDistribution)
        .def("SetDistribution", &DistributedModel<3>::SetDistribution)
        .def("GetTranspose", &DistributedModel<3>::GetTranspose)
        .def("ForceSingleProc", &DistributedModel<3>::ForceSingleProc)
        .def("SetupMPI", &DistributedModel<3>::SetupMPI)
        .def("CreateInitialShape", &DistributedModel<3>::CreateInitialShape)
        .def("GetGlobalShape", &DistributedModel<3>::GetGlobalShape)
        .def("GetLocalStartIndex", DistributedModel_3___GetLocalStartIndexint_int)
        .def("GetLocalIndexRange", DistributedModel_3___GetLocalIndexRangeint_int)
        .def("ApplyConfigSection", &DistributedModel<3>::ApplyConfigSection)
        .def("GetLocalStartIndex", DistributedModel_3___GetLocalStartIndexint_int_int)
        .def("GetLocalIndexRange", DistributedModel_3___GetLocalIndexRangeint_int_int)
        .def("GetGlobalSum", DistributedModel_3___GetGlobalSumdouble)
        .def("GetGlobalSum", DistributedModel_3___GetGlobalSumstd__complex_double_)
        .def("GetGlobalSum", DistributedModel_3___GetGlobalSumblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1__)
        .def("GlobalBarrier", &DistributedModel<3>::GlobalBarrier)
        .def("IsDistributedRank", &DistributedModel<3>::IsDistributedRank)
        .def("ChangeDistribution", &DistributedModel<3>::ChangeDistribution)
        .def("InitMPI", &DistributedModel<3>::InitMPI)
        .def("FinalizeMPI", &DistributedModel<3>::FinalizeMPI)
        .staticmethod("InitMPI")
        .staticmethod("FinalizeMPI")
        .staticmethod("ForceSingleProc")
    );
    register_ptr_to_python< boost::shared_ptr< DistributedModel<3> > >();
    delete DistributedModel_3_scope;

    scope* DistributedModel_4_scope = new scope(
    class_< DistributedModel<4> >("DistributedModel_4", init<  >())
        .def(init< boost::shared_ptr<Distribution> >())
        .def(init< const DistributedModel<4>& >())
        .def_readwrite("ProcId", &DistributedModel<4>::ProcId)
        .def_readwrite("ProcCount", &DistributedModel<4>::ProcCount)
        .def("CreateSubDistributedModel", &DistributedModel<4>::CreateSubDistributedModel)
        .def("IsSingleProc", &DistributedModel<4>::IsSingleProc)
        .def("IsFirstProc", &DistributedModel<4>::IsFirstProc)
        .def("IsLastProc", &DistributedModel<4>::IsLastProc)
        .def("GetDistribution", &DistributedModel<4>::GetDistribution)
        .def("SetDistribution", &DistributedModel<4>::SetDistribution)
        .def("GetTranspose", &DistributedModel<4>::GetTranspose)
        .def("ForceSingleProc", &DistributedModel<4>::ForceSingleProc)
        .def("SetupMPI", &DistributedModel<4>::SetupMPI)
        .def("CreateInitialShape", &DistributedModel<4>::CreateInitialShape)
        .def("GetGlobalShape", &DistributedModel<4>::GetGlobalShape)
        .def("GetLocalStartIndex", DistributedModel_4___GetLocalStartIndexint_int)
        .def("GetLocalIndexRange", DistributedModel_4___GetLocalIndexRangeint_int)
        .def("ApplyConfigSection", &DistributedModel<4>::ApplyConfigSection)
        .def("GetLocalStartIndex", DistributedModel_4___GetLocalStartIndexint_int_int)
        .def("GetLocalIndexRange", DistributedModel_4___GetLocalIndexRangeint_int_int)
        .def("GetGlobalSum", DistributedModel_4___GetGlobalSumdouble)
        .def("GetGlobalSum", DistributedModel_4___GetGlobalSumstd__complex_double_)
        .def("GetGlobalSum", DistributedModel_4___GetGlobalSumblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1__)
        .def("GlobalBarrier", &DistributedModel<4>::GlobalBarrier)
        .def("IsDistributedRank", &DistributedModel<4>::IsDistributedRank)
        .def("ChangeDistribution", &DistributedModel<4>::ChangeDistribution)
        .def("InitMPI", &DistributedModel<4>::InitMPI)
        .def("FinalizeMPI", &DistributedModel<4>::FinalizeMPI)
        .staticmethod("InitMPI")
        .staticmethod("FinalizeMPI")
        .staticmethod("ForceSingleProc")
    );
    register_ptr_to_python< boost::shared_ptr< DistributedModel<4> > >();
    delete DistributedModel_4_scope;

}


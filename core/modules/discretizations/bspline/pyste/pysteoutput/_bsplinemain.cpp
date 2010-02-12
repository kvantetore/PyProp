
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <bspline.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

double (BSpline::BSpline::*BSpline__BSpline__BSplineOverlapIntegralint_int)(int, int)  = &BSpline::BSpline::BSplineOverlapIntegral;

double (BSpline::BSpline::*BSpline__BSpline__BSplineOverlapIntegralblitz__Array_double_1__int_int)(blitz::Array<double,1>, int, int)  = &BSpline::BSpline::BSplineOverlapIntegral;

blitz::Array<std::complex<double>,1> (BSpline::BSpline::*BSpline__BSpline__ExpandFunctionInBSplinesboost__python__object)(boost::python::object)  = &BSpline::BSpline::ExpandFunctionInBSplines;

void (BSpline::BSpline::*BSpline__BSpline__ExpandFunctionInBSplinesblitz__Array_std__complex_double__1__blitz__Array_std__complex_double__1_)(blitz::Array<std::complex<double>,1>, blitz::Array<std::complex<double>,1>)  = &BSpline::BSpline::ExpandFunctionInBSplines;

blitz::Array<std::complex<double>,1> (BSpline::BSpline::*BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1_)(blitz::Array<std::complex<double>,1>)  = &BSpline::BSpline::ConstructFunctionFromBSplineExpansion;

blitz::Array<std::complex<double>,1> (BSpline::BSpline::*BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1__blitz__Array_double_1_)(blitz::Array<std::complex<double>,1>, blitz::Array<double,1>)  = &BSpline::BSpline::ConstructFunctionFromBSplineExpansion;

void (BSpline::BSpline::*BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1__blitz__Array_double_1__blitz__Array_std__complex_double__1_)(blitz::Array<std::complex<double>,1>, blitz::Array<double,1>, blitz::Array<std::complex<double>,1>)  = &BSpline::BSpline::ConstructFunctionFromBSplineExpansion;


}// namespace 


// Module ======================================================================
void Export_core_modules_discretizations_bspline_pyste_bsplinemain()
{
    scope* BSpline_BSpline_scope = new scope(
    class_< BSpline::BSpline, boost::noncopyable >("BSpline", init<  >())
        .def_readwrite("NumberOfBSplines", &BSpline::BSpline::NumberOfBSplines)
        .def_readwrite("MaxSplineOrder", &BSpline::BSpline::MaxSplineOrder)
        .def_readwrite("ProjectionAlgorithm", &BSpline::BSpline::ProjectionAlgorithm)
        .def_readwrite("LapackAlgorithm", &BSpline::BSpline::LapackAlgorithm)
        .def("GetBreakpointSequence", &BSpline::BSpline::GetBreakpointSequence)
        .def("GetContinuitySequence", &BSpline::BSpline::GetContinuitySequence)
        .def("GetKnotSequence", &BSpline::BSpline::GetKnotSequence)
        .def("GetKnotGridIndexMap", &BSpline::BSpline::GetKnotGridIndexMap)
        .def("GetTopKnotMap", &BSpline::BSpline::GetTopKnotMap)
        .def("GetWeights", &BSpline::BSpline::GetWeights)
        .def("GetNodes", &BSpline::BSpline::GetNodes)
        .def("GetQuadratureGrid", &BSpline::BSpline::GetQuadratureGrid)
        .def("GetQuadratureGridGlobal", &BSpline::BSpline::GetQuadratureGridGlobal)
        .def("GetQuadratureWeightsGlobal", &BSpline::BSpline::GetQuadratureWeightsGlobal)
        .def("GetOverlap", &BSpline::BSpline::GetOverlap)
        .def("ResizeBreakpointSequence", &BSpline::BSpline::ResizeBreakpointSequence)
        .def("ResizeContinuitySequence", &BSpline::BSpline::ResizeContinuitySequence)
        .def("ResizeKnotSequence", &BSpline::BSpline::ResizeKnotSequence)
        .def("ResizeKnotGridIndexMap", &BSpline::BSpline::ResizeKnotGridIndexMap)
        .def("ResizeTopKnotMap", &BSpline::BSpline::ResizeTopKnotMap)
        .def("ResizeWeights", &BSpline::BSpline::ResizeWeights)
        .def("ResizeNodes", &BSpline::BSpline::ResizeNodes)
        .def("ResizeQuadratureGridGlobal", &BSpline::BSpline::ResizeQuadratureGridGlobal)
        .def("ResizeQuadratureWeightsGlobal", &BSpline::BSpline::ResizeQuadratureWeightsGlobal)
        .def("EvaluateBSpline", &BSpline::BSpline::EvaluateBSpline)
        .def("EvaluateBSplineDerivative1", &BSpline::BSpline::EvaluateBSplineDerivative1)
        .def("EvaluateBSplineDerivative2", &BSpline::BSpline::EvaluateBSplineDerivative2)
        .def("EvaluateBSplineOnGrid", &BSpline::BSpline::EvaluateBSplineOnGrid)
        .def("GetBSpline", &BSpline::BSpline::GetBSpline)
        .def("GetBSplineDerivative1", &BSpline::BSpline::GetBSplineDerivative1)
        .def("GetBSplineDerivative2", &BSpline::BSpline::GetBSplineDerivative2)
        .def("BSplineOverlapIntegral", BSpline__BSpline__BSplineOverlapIntegralint_int)
        .def("BSplineOverlapIntegral", BSpline__BSpline__BSplineOverlapIntegralblitz__Array_double_1__int_int)
        .def("ProjectOnBSpline", &BSpline::BSpline::ProjectOnBSpline)
        .def("BSplineDerivative2OverlapIntegral", &BSpline::BSpline::BSplineDerivative2OverlapIntegral)
        .def("SolveForOverlapMatrix", &BSpline::BSpline::SolveForOverlapMatrix)
        .def("BSplineGlobalOverlapIntegral", &BSpline::BSpline::BSplineGlobalOverlapIntegral_double)
        .def("BSplineGlobalOverlapIntegral", &BSpline::BSpline::BSplineGlobalOverlapIntegral_cplx)
        .def("CreateBSplineTable", &BSpline::BSpline::CreateBSplineTable)
        .def("CreateBSplineDerivativeTable", &BSpline::BSpline::CreateBSplineDerivativeTable)
        .def("SetupOverlap", &BSpline::BSpline::SetupOverlap)
        .def("ExpandFunctionInBSplines", BSpline__BSpline__ExpandFunctionInBSplinesboost__python__object)
        .def("ExpandFunctionInBSplines", BSpline__BSpline__ExpandFunctionInBSplinesblitz__Array_std__complex_double__1__blitz__Array_std__complex_double__1_)
        .def("ConstructFunctionFromBSplineExpansion", BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1_)
        .def("ConstructFunctionFromBSplineExpansion", BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1__blitz__Array_double_1_)
        .def("ConstructFunctionFromBSplineExpansion", BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1__blitz__Array_double_1__blitz__Array_std__complex_double__1_)
        .def("ComputeStartIndex", &BSpline::BSpline::ComputeStartIndex)
        .def("ComputeStopIndex", &BSpline::BSpline::ComputeStopIndex)
        .def("ComputeIndexShift", &BSpline::BSpline::ComputeIndexShift)
        .def("GetGridIndex", &BSpline::BSpline::GetGridIndex)
        .def("ScaleAndTranslate", &BSpline::BSpline::ScaleAndTranslate)
    );
    register_ptr_to_python< boost::shared_ptr< BSpline::BSpline > >();
    delete BSpline_BSpline_scope;

}


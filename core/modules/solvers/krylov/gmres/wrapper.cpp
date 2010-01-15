
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <gmreswrapper.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libgmres)
{
    class_< krylov::GmresWrapper<1>, boost::noncopyable >("krylov_GmresWrapper_1", init<  >())
        .def("ApplyConfigSection", &krylov::GmresWrapper<1>::ApplyConfigSection)
        .def("Setup", &krylov::GmresWrapper<1>::Setup)
        .def("Solve", &krylov::GmresWrapper<1>::Solve)
        .def("ApplyOperator", &krylov::GmresWrapper<1>::ApplyOperator)
        .def("EstimateMemoryUsage", &krylov::GmresWrapper<1>::EstimateMemoryUsage)
        .def("GetOperatorCount", &krylov::GmresWrapper<1>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::GmresWrapper<1>::GetOrthogonalizationCount)
        .def("PrintStatistics", &krylov::GmresWrapper<1>::PrintStatistics)
        .def("ResetStatistics", &krylov::GmresWrapper<1>::ResetStatistics)
        .def("GetErrorEstimate", &krylov::GmresWrapper<1>::GetErrorEstimate)
        .def("GetErrorEstimateList", &krylov::GmresWrapper<1>::GetErrorEstimateList)
        .def("GetHessenbergMatrix", &krylov::GmresWrapper<1>::GetHessenbergMatrix)
        .def("SetupResidual", &PypropKrylovWrapper::SetupResidual)
    ;

    class_< krylov::GmresWrapper<2>, boost::noncopyable >("krylov_GmresWrapper_2", init<  >())
        .def("ApplyConfigSection", &krylov::GmresWrapper<2>::ApplyConfigSection)
        .def("Setup", &krylov::GmresWrapper<2>::Setup)
        .def("Solve", &krylov::GmresWrapper<2>::Solve)
        .def("ApplyOperator", &krylov::GmresWrapper<2>::ApplyOperator)
        .def("EstimateMemoryUsage", &krylov::GmresWrapper<2>::EstimateMemoryUsage)
        .def("GetOperatorCount", &krylov::GmresWrapper<2>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::GmresWrapper<2>::GetOrthogonalizationCount)
        .def("PrintStatistics", &krylov::GmresWrapper<2>::PrintStatistics)
        .def("ResetStatistics", &krylov::GmresWrapper<2>::ResetStatistics)
        .def("GetErrorEstimate", &krylov::GmresWrapper<2>::GetErrorEstimate)
        .def("GetErrorEstimateList", &krylov::GmresWrapper<2>::GetErrorEstimateList)
        .def("GetHessenbergMatrix", &krylov::GmresWrapper<2>::GetHessenbergMatrix)
        .def("SetupResidual", &PypropKrylovWrapper::SetupResidual)
    ;

    class_< krylov::GmresWrapper<3>, boost::noncopyable >("krylov_GmresWrapper_3", init<  >())
        .def("ApplyConfigSection", &krylov::GmresWrapper<3>::ApplyConfigSection)
        .def("Setup", &krylov::GmresWrapper<3>::Setup)
        .def("Solve", &krylov::GmresWrapper<3>::Solve)
        .def("ApplyOperator", &krylov::GmresWrapper<3>::ApplyOperator)
        .def("EstimateMemoryUsage", &krylov::GmresWrapper<3>::EstimateMemoryUsage)
        .def("GetOperatorCount", &krylov::GmresWrapper<3>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::GmresWrapper<3>::GetOrthogonalizationCount)
        .def("PrintStatistics", &krylov::GmresWrapper<3>::PrintStatistics)
        .def("ResetStatistics", &krylov::GmresWrapper<3>::ResetStatistics)
        .def("GetErrorEstimate", &krylov::GmresWrapper<3>::GetErrorEstimate)
        .def("GetErrorEstimateList", &krylov::GmresWrapper<3>::GetErrorEstimateList)
        .def("GetHessenbergMatrix", &krylov::GmresWrapper<3>::GetHessenbergMatrix)
        .def("SetupResidual", &PypropKrylovWrapper::SetupResidual)
    ;

    class_< krylov::GmresWrapper<4>, boost::noncopyable >("krylov_GmresWrapper_4", init<  >())
        .def("ApplyConfigSection", &krylov::GmresWrapper<4>::ApplyConfigSection)
        .def("Setup", &krylov::GmresWrapper<4>::Setup)
        .def("Solve", &krylov::GmresWrapper<4>::Solve)
        .def("ApplyOperator", &krylov::GmresWrapper<4>::ApplyOperator)
        .def("EstimateMemoryUsage", &krylov::GmresWrapper<4>::EstimateMemoryUsage)
        .def("GetOperatorCount", &krylov::GmresWrapper<4>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::GmresWrapper<4>::GetOrthogonalizationCount)
        .def("PrintStatistics", &krylov::GmresWrapper<4>::PrintStatistics)
        .def("ResetStatistics", &krylov::GmresWrapper<4>::ResetStatistics)
        .def("GetErrorEstimate", &krylov::GmresWrapper<4>::GetErrorEstimate)
        .def("GetErrorEstimateList", &krylov::GmresWrapper<4>::GetErrorEstimateList)
        .def("GetHessenbergMatrix", &krylov::GmresWrapper<4>::GetHessenbergMatrix)
        .def("SetupResidual", &PypropKrylovWrapper::SetupResidual)
    ;

}


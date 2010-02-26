
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <pampwrapper.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpamp)
{
    class_< krylov::PampWrapper<1>, boost::noncopyable >("krylov_PampWrapper_1", init<  >())
        .def("ApplyConfigSection", &krylov::PampWrapper<1>::ApplyConfigSection)
        .def("Setup", &krylov::PampWrapper<1>::Setup)
        .def("AdvanceStep", &krylov::PampWrapper<1>::AdvanceStep)
        .def("ApplyOperator", &krylov::PampWrapper<1>::ApplyOperator)
        .def("EstimateMemoryUsage", &krylov::PampWrapper<1>::EstimateMemoryUsage)
        .def("GetOperatorCount", &krylov::PampWrapper<1>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::PampWrapper<1>::GetOrthogonalizationCount)
        .def("PrintStatistics", &krylov::PampWrapper<1>::PrintStatistics)
        .def("ResetStatistics", &krylov::PampWrapper<1>::ResetStatistics)
        .def("TestPade", &krylov::PampWrapper<1>::TestPade)
        .def("GetResidualNorm", &krylov::PampWrapper<1>::GetResidualNorm)
        .def("GetHessenbergMatrix", &krylov::PampWrapper<1>::GetHessenbergMatrix)
        .def("GetHessenbergMatrixExp", &krylov::PampWrapper<1>::GetHessenbergMatrixExp)
        .def("GetPropagationErrorEstimate", &krylov::PampWrapper<1>::GetPropagationErrorEstimate)
        .def("SetupResidual", &PypropKrylovWrapper::SetupResidual)
    ;

    class_< krylov::PampWrapper<2>, boost::noncopyable >("krylov_PampWrapper_2", init<  >())
        .def("ApplyConfigSection", &krylov::PampWrapper<2>::ApplyConfigSection)
        .def("Setup", &krylov::PampWrapper<2>::Setup)
        .def("AdvanceStep", &krylov::PampWrapper<2>::AdvanceStep)
        .def("ApplyOperator", &krylov::PampWrapper<2>::ApplyOperator)
        .def("EstimateMemoryUsage", &krylov::PampWrapper<2>::EstimateMemoryUsage)
        .def("GetOperatorCount", &krylov::PampWrapper<2>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::PampWrapper<2>::GetOrthogonalizationCount)
        .def("PrintStatistics", &krylov::PampWrapper<2>::PrintStatistics)
        .def("ResetStatistics", &krylov::PampWrapper<2>::ResetStatistics)
        .def("TestPade", &krylov::PampWrapper<2>::TestPade)
        .def("GetResidualNorm", &krylov::PampWrapper<2>::GetResidualNorm)
        .def("GetHessenbergMatrix", &krylov::PampWrapper<2>::GetHessenbergMatrix)
        .def("GetHessenbergMatrixExp", &krylov::PampWrapper<2>::GetHessenbergMatrixExp)
        .def("GetPropagationErrorEstimate", &krylov::PampWrapper<2>::GetPropagationErrorEstimate)
        .def("SetupResidual", &PypropKrylovWrapper::SetupResidual)
    ;

    class_< krylov::PampWrapper<3>, boost::noncopyable >("krylov_PampWrapper_3", init<  >())
        .def("ApplyConfigSection", &krylov::PampWrapper<3>::ApplyConfigSection)
        .def("Setup", &krylov::PampWrapper<3>::Setup)
        .def("AdvanceStep", &krylov::PampWrapper<3>::AdvanceStep)
        .def("ApplyOperator", &krylov::PampWrapper<3>::ApplyOperator)
        .def("EstimateMemoryUsage", &krylov::PampWrapper<3>::EstimateMemoryUsage)
        .def("GetOperatorCount", &krylov::PampWrapper<3>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::PampWrapper<3>::GetOrthogonalizationCount)
        .def("PrintStatistics", &krylov::PampWrapper<3>::PrintStatistics)
        .def("ResetStatistics", &krylov::PampWrapper<3>::ResetStatistics)
        .def("TestPade", &krylov::PampWrapper<3>::TestPade)
        .def("GetResidualNorm", &krylov::PampWrapper<3>::GetResidualNorm)
        .def("GetHessenbergMatrix", &krylov::PampWrapper<3>::GetHessenbergMatrix)
        .def("GetHessenbergMatrixExp", &krylov::PampWrapper<3>::GetHessenbergMatrixExp)
        .def("GetPropagationErrorEstimate", &krylov::PampWrapper<3>::GetPropagationErrorEstimate)
        .def("SetupResidual", &PypropKrylovWrapper::SetupResidual)
    ;

    class_< krylov::PampWrapper<4>, boost::noncopyable >("krylov_PampWrapper_4", init<  >())
        .def("ApplyConfigSection", &krylov::PampWrapper<4>::ApplyConfigSection)
        .def("Setup", &krylov::PampWrapper<4>::Setup)
        .def("AdvanceStep", &krylov::PampWrapper<4>::AdvanceStep)
        .def("ApplyOperator", &krylov::PampWrapper<4>::ApplyOperator)
        .def("EstimateMemoryUsage", &krylov::PampWrapper<4>::EstimateMemoryUsage)
        .def("GetOperatorCount", &krylov::PampWrapper<4>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::PampWrapper<4>::GetOrthogonalizationCount)
        .def("PrintStatistics", &krylov::PampWrapper<4>::PrintStatistics)
        .def("ResetStatistics", &krylov::PampWrapper<4>::ResetStatistics)
        .def("TestPade", &krylov::PampWrapper<4>::TestPade)
        .def("GetResidualNorm", &krylov::PampWrapper<4>::GetResidualNorm)
        .def("GetHessenbergMatrix", &krylov::PampWrapper<4>::GetHessenbergMatrix)
        .def("GetHessenbergMatrixExp", &krylov::PampWrapper<4>::GetHessenbergMatrixExp)
        .def("GetPropagationErrorEstimate", &krylov::PampWrapper<4>::GetPropagationErrorEstimate)
        .def("SetupResidual", &PypropKrylovWrapper::SetupResidual)
    ;

}


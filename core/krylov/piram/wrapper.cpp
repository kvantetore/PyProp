
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <piramsolver.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpiram)
{
    class_< krylov::PiramSolver<1>, boost::noncopyable >("krylov_PiramSolver_1", init<  >())
        .def("ApplyConfigSection", &krylov::PiramSolver<1>::ApplyConfigSection)
        .def("Setup", &krylov::PiramSolver<1>::Setup)
        .def("Solve", &krylov::PiramSolver<1>::Solve)
        .def("SetupResidual", &krylov::PiramSolver<1>::SetupResidual)
        .def("ApplyOperator", &krylov::PiramSolver<1>::ApplyOperator)
        .def("GetEigenvalues", &krylov::PiramSolver<1>::GetEigenvalues)
        .def("GetEigenvector", &krylov::PiramSolver<1>::GetEigenvector)
        .def("EstimateMemoryUsage", &krylov::PiramSolver<1>::EstimateMemoryUsage)
        .def("GetErrorEstimates", &krylov::PiramSolver<1>::GetErrorEstimates)
        .def("GetConvergenceEstimates", &krylov::PiramSolver<1>::GetConvergenceEstimates)
        .def("GetEigenvalueCount", &krylov::PiramSolver<1>::GetEigenvalueCount)
        .def("GetRestartCount", &krylov::PiramSolver<1>::GetRestartCount)
        .def("GetOperatorCount", &krylov::PiramSolver<1>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::PiramSolver<1>::GetOrthogonalizationCount)
    ;

    class_< krylov::PiramSolver<2>, boost::noncopyable >("krylov_PiramSolver_2", init<  >())
        .def("ApplyConfigSection", &krylov::PiramSolver<2>::ApplyConfigSection)
        .def("Setup", &krylov::PiramSolver<2>::Setup)
        .def("Solve", &krylov::PiramSolver<2>::Solve)
        .def("SetupResidual", &krylov::PiramSolver<2>::SetupResidual)
        .def("ApplyOperator", &krylov::PiramSolver<2>::ApplyOperator)
        .def("GetEigenvalues", &krylov::PiramSolver<2>::GetEigenvalues)
        .def("GetEigenvector", &krylov::PiramSolver<2>::GetEigenvector)
        .def("EstimateMemoryUsage", &krylov::PiramSolver<2>::EstimateMemoryUsage)
        .def("GetErrorEstimates", &krylov::PiramSolver<2>::GetErrorEstimates)
        .def("GetConvergenceEstimates", &krylov::PiramSolver<2>::GetConvergenceEstimates)
        .def("GetEigenvalueCount", &krylov::PiramSolver<2>::GetEigenvalueCount)
        .def("GetRestartCount", &krylov::PiramSolver<2>::GetRestartCount)
        .def("GetOperatorCount", &krylov::PiramSolver<2>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::PiramSolver<2>::GetOrthogonalizationCount)
    ;

    class_< krylov::PiramSolver<3>, boost::noncopyable >("krylov_PiramSolver_3", init<  >())
        .def("ApplyConfigSection", &krylov::PiramSolver<3>::ApplyConfigSection)
        .def("Setup", &krylov::PiramSolver<3>::Setup)
        .def("Solve", &krylov::PiramSolver<3>::Solve)
        .def("SetupResidual", &krylov::PiramSolver<3>::SetupResidual)
        .def("ApplyOperator", &krylov::PiramSolver<3>::ApplyOperator)
        .def("GetEigenvalues", &krylov::PiramSolver<3>::GetEigenvalues)
        .def("GetEigenvector", &krylov::PiramSolver<3>::GetEigenvector)
        .def("EstimateMemoryUsage", &krylov::PiramSolver<3>::EstimateMemoryUsage)
        .def("GetErrorEstimates", &krylov::PiramSolver<3>::GetErrorEstimates)
        .def("GetConvergenceEstimates", &krylov::PiramSolver<3>::GetConvergenceEstimates)
        .def("GetEigenvalueCount", &krylov::PiramSolver<3>::GetEigenvalueCount)
        .def("GetRestartCount", &krylov::PiramSolver<3>::GetRestartCount)
        .def("GetOperatorCount", &krylov::PiramSolver<3>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::PiramSolver<3>::GetOrthogonalizationCount)
    ;

    class_< krylov::PiramSolver<4>, boost::noncopyable >("krylov_PiramSolver_4", init<  >())
        .def("ApplyConfigSection", &krylov::PiramSolver<4>::ApplyConfigSection)
        .def("Setup", &krylov::PiramSolver<4>::Setup)
        .def("Solve", &krylov::PiramSolver<4>::Solve)
        .def("SetupResidual", &krylov::PiramSolver<4>::SetupResidual)
        .def("ApplyOperator", &krylov::PiramSolver<4>::ApplyOperator)
        .def("GetEigenvalues", &krylov::PiramSolver<4>::GetEigenvalues)
        .def("GetEigenvector", &krylov::PiramSolver<4>::GetEigenvector)
        .def("EstimateMemoryUsage", &krylov::PiramSolver<4>::EstimateMemoryUsage)
        .def("GetErrorEstimates", &krylov::PiramSolver<4>::GetErrorEstimates)
        .def("GetConvergenceEstimates", &krylov::PiramSolver<4>::GetConvergenceEstimates)
        .def("GetEigenvalueCount", &krylov::PiramSolver<4>::GetEigenvalueCount)
        .def("GetRestartCount", &krylov::PiramSolver<4>::GetRestartCount)
        .def("GetOperatorCount", &krylov::PiramSolver<4>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &krylov::PiramSolver<4>::GetOrthogonalizationCount)
    ;

}


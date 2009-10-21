
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <anasaziwrapper.h>
#include <ifpackwrapper.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct AnasaziSolver_1_Wrapper: AnasaziSolver<1>
{
    AnasaziSolver_1_Wrapper(PyObject* py_self_):
        AnasaziSolver<1>(), py_self(py_self_) {}

    void SetupResidual(blitz::Array<std::complex<double>,1>& p0) {
        call_method< void >(py_self, "SetupResidual", p0);
    }

    void default_SetupResidual(blitz::Array<std::complex<double>,1>& p0) {
        AnasaziSolver<1>::SetupResidual(p0);
    }

    PyObject* py_self;
};

struct AnasaziSolver_2_Wrapper: AnasaziSolver<2>
{
    AnasaziSolver_2_Wrapper(PyObject* py_self_):
        AnasaziSolver<2>(), py_self(py_self_) {}

    void SetupResidual(blitz::Array<std::complex<double>,1>& p0) {
        call_method< void >(py_self, "SetupResidual", p0);
    }

    void default_SetupResidual(blitz::Array<std::complex<double>,1>& p0) {
        AnasaziSolver<2>::SetupResidual(p0);
    }

    PyObject* py_self;
};

struct AnasaziSolver_3_Wrapper: AnasaziSolver<3>
{
    AnasaziSolver_3_Wrapper(PyObject* py_self_):
        AnasaziSolver<3>(), py_self(py_self_) {}

    void SetupResidual(blitz::Array<std::complex<double>,1>& p0) {
        call_method< void >(py_self, "SetupResidual", p0);
    }

    void default_SetupResidual(blitz::Array<std::complex<double>,1>& p0) {
        AnasaziSolver<3>::SetupResidual(p0);
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libtrilinos)
{
    class_< AnasaziSolver<1>, boost::noncopyable, AnasaziSolver_1_Wrapper >("AnasaziSolver_1", init<  >())
        .def("SetupResidual", &AnasaziSolver<1>::SetupResidual, &AnasaziSolver_1_Wrapper::default_SetupResidual)
        .def("ApplyConfigSection", &AnasaziSolver<1>::ApplyConfigSection)
        .def("Setup", &AnasaziSolver<1>::Setup)
        .def("Solve", &AnasaziSolver<1>::Solve)
        .def("GetEigenvalues", &AnasaziSolver<1>::GetEigenvalues)
        .def("GetEigenvector", &AnasaziSolver<1>::GetEigenvector)
        .def("EstimateMemoryUsage", &AnasaziSolver<1>::EstimateMemoryUsage)
        .def("GetErrorEstimates", &AnasaziSolver<1>::GetErrorEstimates)
        .def("GetConvergenceEstimates", &AnasaziSolver<1>::GetConvergenceEstimates)
        .def("GetEigenvalueCount", &AnasaziSolver<1>::GetEigenvalueCount)
        .def("GetRestartCount", &AnasaziSolver<1>::GetRestartCount)
        .def("GetOperatorCount", &AnasaziSolver<1>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &AnasaziSolver<1>::GetOrthogonalizationCount)
    ;

    class_< AnasaziSolver<2>, boost::noncopyable, AnasaziSolver_2_Wrapper >("AnasaziSolver_2", init<  >())
        .def("SetupResidual", &AnasaziSolver<2>::SetupResidual, &AnasaziSolver_2_Wrapper::default_SetupResidual)
        .def("ApplyConfigSection", &AnasaziSolver<2>::ApplyConfigSection)
        .def("Setup", &AnasaziSolver<2>::Setup)
        .def("Solve", &AnasaziSolver<2>::Solve)
        .def("GetEigenvalues", &AnasaziSolver<2>::GetEigenvalues)
        .def("GetEigenvector", &AnasaziSolver<2>::GetEigenvector)
        .def("EstimateMemoryUsage", &AnasaziSolver<2>::EstimateMemoryUsage)
        .def("GetErrorEstimates", &AnasaziSolver<2>::GetErrorEstimates)
        .def("GetConvergenceEstimates", &AnasaziSolver<2>::GetConvergenceEstimates)
        .def("GetEigenvalueCount", &AnasaziSolver<2>::GetEigenvalueCount)
        .def("GetRestartCount", &AnasaziSolver<2>::GetRestartCount)
        .def("GetOperatorCount", &AnasaziSolver<2>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &AnasaziSolver<2>::GetOrthogonalizationCount)
    ;

    class_< AnasaziSolver<3>, boost::noncopyable, AnasaziSolver_3_Wrapper >("AnasaziSolver_3", init<  >())
        .def("SetupResidual", &AnasaziSolver<3>::SetupResidual, &AnasaziSolver_3_Wrapper::default_SetupResidual)
        .def("ApplyConfigSection", &AnasaziSolver<3>::ApplyConfigSection)
        .def("Setup", &AnasaziSolver<3>::Setup)
        .def("Solve", &AnasaziSolver<3>::Solve)
        .def("GetEigenvalues", &AnasaziSolver<3>::GetEigenvalues)
        .def("GetEigenvector", &AnasaziSolver<3>::GetEigenvector)
        .def("EstimateMemoryUsage", &AnasaziSolver<3>::EstimateMemoryUsage)
        .def("GetErrorEstimates", &AnasaziSolver<3>::GetErrorEstimates)
        .def("GetConvergenceEstimates", &AnasaziSolver<3>::GetConvergenceEstimates)
        .def("GetEigenvalueCount", &AnasaziSolver<3>::GetEigenvalueCount)
        .def("GetRestartCount", &AnasaziSolver<3>::GetRestartCount)
        .def("GetOperatorCount", &AnasaziSolver<3>::GetOperatorCount)
        .def("GetOrthogonalizationCount", &AnasaziSolver<3>::GetOrthogonalizationCount)
    ;

    class_< IfpackRadialPreconditioner<1>, boost::noncopyable >("IfpackRadialPreconditioner_1", init<  >())
        .def("ApplyConfigSection", &IfpackRadialPreconditioner<1>::ApplyConfigSection)
        .def("Setup", &IfpackRadialPreconditioner<1>::Setup)
        .def("Solve", &IfpackRadialPreconditioner<1>::Solve)
    ;

    class_< IfpackRadialPreconditioner<2>, boost::noncopyable >("IfpackRadialPreconditioner_2", init<  >())
        .def("ApplyConfigSection", &IfpackRadialPreconditioner<2>::ApplyConfigSection)
        .def("Setup", &IfpackRadialPreconditioner<2>::Setup)
        .def("Solve", &IfpackRadialPreconditioner<2>::Solve)
    ;

    class_< IfpackRadialPreconditioner<3>, boost::noncopyable >("IfpackRadialPreconditioner_3", init<  >())
        .def("ApplyConfigSection", &IfpackRadialPreconditioner<3>::ApplyConfigSection)
        .def("Setup", &IfpackRadialPreconditioner<3>::Setup)
        .def("Solve", &IfpackRadialPreconditioner<3>::Solve)
    ;

}


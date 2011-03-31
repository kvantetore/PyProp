
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <finitediff/cranknicholsonpropagator.h>
#include <finitediff/finitedifferencehelper.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct FiniteDifferenceHelper_Wrapper: FiniteDifferenceHelper
{
    FiniteDifferenceHelper_Wrapper(PyObject* py_self_, const FiniteDifferenceHelper& p0):
        FiniteDifferenceHelper(p0), py_self(py_self_) {}

    FiniteDifferenceHelper_Wrapper(PyObject* py_self_):
        FiniteDifferenceHelper(), py_self(py_self_) {}

    blitz::Array<std::complex<double>,1> FindDifferenceCoefficients(int p0) {
        return call_method< blitz::Array<std::complex<double>,1> >(py_self, "FindDifferenceCoefficients", p0);
    }

    blitz::Array<std::complex<double>,1> default_FindDifferenceCoefficients(int p0) {
        return FiniteDifferenceHelper::FindDifferenceCoefficients(p0);
    }

    PyObject* py_self;
};

struct FiniteDifferenceHelperCustomBoundary_Wrapper: FiniteDifferenceHelperCustomBoundary
{
    FiniteDifferenceHelperCustomBoundary_Wrapper(PyObject* py_self_, const FiniteDifferenceHelperCustomBoundary& p0):
        FiniteDifferenceHelperCustomBoundary(p0), py_self(py_self_) {}

    FiniteDifferenceHelperCustomBoundary_Wrapper(PyObject* py_self_):
        FiniteDifferenceHelperCustomBoundary(), py_self(py_self_) {}

    blitz::Array<std::complex<double>,1> FindDifferenceCoefficients(int p0) {
        return call_method< blitz::Array<std::complex<double>,1> >(py_self, "FindDifferenceCoefficients", p0);
    }

    blitz::Array<std::complex<double>,1> default_FindDifferenceCoefficients(int p0) {
        return FiniteDifferenceHelperCustomBoundary::FindDifferenceCoefficients(p0);
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
void Export_python_cranknicholson()
{
    class_< CrankNicholsonPropagator<1> >("CrankNicholsonPropagator_1", init<  >())
        .def(init< const CrankNicholsonPropagator<1>& >())
        .def("MultiplyKineticEnergyOperator", &CrankNicholsonPropagator<1>::MultiplyKineticEnergyOperator)
        .def("AdvanceStep", &CrankNicholsonPropagator<1>::AdvanceStep)
        .def("SetupStep", &CrankNicholsonPropagator<1>::SetupStep)
        .def("ApplyConfigSection", &CrankNicholsonPropagator<1>::ApplyConfigSection)
        .def("GetLaplacianBlasBanded", &CrankNicholsonPropagator<1>::GetLaplacianBlasBanded)
        .def("GetLaplacianDistributedBanded", &CrankNicholsonPropagator<1>::GetLaplacianDistributedBanded)
        .def("GetLaplacianFull", &CrankNicholsonPropagator<1>::GetLaplacianFull)
        .def("GetBackwardPropagationLapackBanded", &CrankNicholsonPropagator<1>::GetBackwardPropagationLapackBanded)
    ;

    class_< CrankNicholsonPropagator<2> >("CrankNicholsonPropagator_2", init<  >())
        .def(init< const CrankNicholsonPropagator<2>& >())
        .def("MultiplyKineticEnergyOperator", &CrankNicholsonPropagator<2>::MultiplyKineticEnergyOperator)
        .def("AdvanceStep", &CrankNicholsonPropagator<2>::AdvanceStep)
        .def("SetupStep", &CrankNicholsonPropagator<2>::SetupStep)
        .def("ApplyConfigSection", &CrankNicholsonPropagator<2>::ApplyConfigSection)
        .def("GetLaplacianBlasBanded", &CrankNicholsonPropagator<2>::GetLaplacianBlasBanded)
        .def("GetLaplacianDistributedBanded", &CrankNicholsonPropagator<2>::GetLaplacianDistributedBanded)
        .def("GetLaplacianFull", &CrankNicholsonPropagator<2>::GetLaplacianFull)
        .def("GetBackwardPropagationLapackBanded", &CrankNicholsonPropagator<2>::GetBackwardPropagationLapackBanded)
    ;

    class_< CrankNicholsonPropagator<3> >("CrankNicholsonPropagator_3", init<  >())
        .def(init< const CrankNicholsonPropagator<3>& >())
        .def("MultiplyKineticEnergyOperator", &CrankNicholsonPropagator<3>::MultiplyKineticEnergyOperator)
        .def("AdvanceStep", &CrankNicholsonPropagator<3>::AdvanceStep)
        .def("SetupStep", &CrankNicholsonPropagator<3>::SetupStep)
        .def("ApplyConfigSection", &CrankNicholsonPropagator<3>::ApplyConfigSection)
        .def("GetLaplacianBlasBanded", &CrankNicholsonPropagator<3>::GetLaplacianBlasBanded)
        .def("GetLaplacianDistributedBanded", &CrankNicholsonPropagator<3>::GetLaplacianDistributedBanded)
        .def("GetLaplacianFull", &CrankNicholsonPropagator<3>::GetLaplacianFull)
        .def("GetBackwardPropagationLapackBanded", &CrankNicholsonPropagator<3>::GetBackwardPropagationLapackBanded)
    ;

    class_< CrankNicholsonPropagator<4> >("CrankNicholsonPropagator_4", init<  >())
        .def(init< const CrankNicholsonPropagator<4>& >())
        .def("MultiplyKineticEnergyOperator", &CrankNicholsonPropagator<4>::MultiplyKineticEnergyOperator)
        .def("AdvanceStep", &CrankNicholsonPropagator<4>::AdvanceStep)
        .def("SetupStep", &CrankNicholsonPropagator<4>::SetupStep)
        .def("ApplyConfigSection", &CrankNicholsonPropagator<4>::ApplyConfigSection)
        .def("GetLaplacianBlasBanded", &CrankNicholsonPropagator<4>::GetLaplacianBlasBanded)
        .def("GetLaplacianDistributedBanded", &CrankNicholsonPropagator<4>::GetLaplacianDistributedBanded)
        .def("GetLaplacianFull", &CrankNicholsonPropagator<4>::GetLaplacianFull)
        .def("GetBackwardPropagationLapackBanded", &CrankNicholsonPropagator<4>::GetBackwardPropagationLapackBanded)
    ;

    class_< FiniteDifferenceHelper, FiniteDifferenceHelper_Wrapper >("FiniteDifferenceHelper", init<  >())
        .def(init< const FiniteDifferenceHelper& >())
        .def("FindDifferenceCoefficients", &FiniteDifferenceHelper::FindDifferenceCoefficients, &FiniteDifferenceHelper_Wrapper::default_FindDifferenceCoefficients)
        .def("Setup", &FiniteDifferenceHelper::Setup)
        .def("GetDifferenceOrder", &FiniteDifferenceHelper::GetDifferenceOrder)
        .def("SetupLaplacianBlasBanded", &FiniteDifferenceHelper::SetupLaplacianBlasBanded)
    ;

    class_< FiniteDifferenceHelperCustomBoundary, bases< FiniteDifferenceHelper > , FiniteDifferenceHelperCustomBoundary_Wrapper >("FiniteDifferenceHelperCustomBoundary", init<  >())
        .def(init< const FiniteDifferenceHelperCustomBoundary& >())
        .def("FindDifferenceCoefficients", (blitz::Array<std::complex<double>,1> (FiniteDifferenceHelperCustomBoundary::*)(int) )&FiniteDifferenceHelperCustomBoundary::FindDifferenceCoefficients, (blitz::Array<std::complex<double>,1> (FiniteDifferenceHelperCustomBoundary_Wrapper::*)(int))&FiniteDifferenceHelperCustomBoundary_Wrapper::default_FindDifferenceCoefficients)
        .def("Setup", &FiniteDifferenceHelperCustomBoundary::Setup)
    ;

}


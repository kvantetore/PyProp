
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <src/finitedifferencehelper.h>

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
void Export_pyprop_modules_discretizations_finitedifference_src_pyste_finitedifferencehelper()
{
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


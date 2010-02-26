
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <odewrapper.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libode)
{
    class_< ODE::OdeWrapper<1>, boost::noncopyable >("ODE_OdeWrapper_1", init<  >())
        .def_readwrite("ImTime", &ODE::OdeWrapper<1>::ImTime)
        .def("ApplyConfigSection", &ODE::OdeWrapper<1>::ApplyConfigSection)
        .def("Setup", &ODE::OdeWrapper<1>::Setup)
        .def("AdvanceStep", &ODE::OdeWrapper<1>::AdvanceStep)
        .def("GetPropagatedTime", &ODE::OdeWrapper<1>::GetPropagatedTime)
        .def("SetStartTime", &ODE::OdeWrapper<1>::SetStartTime)
    ;

    class_< ODE::OdeWrapper<2>, boost::noncopyable >("ODE_OdeWrapper_2", init<  >())
        .def_readwrite("ImTime", &ODE::OdeWrapper<2>::ImTime)
        .def("ApplyConfigSection", &ODE::OdeWrapper<2>::ApplyConfigSection)
        .def("Setup", &ODE::OdeWrapper<2>::Setup)
        .def("AdvanceStep", &ODE::OdeWrapper<2>::AdvanceStep)
        .def("GetPropagatedTime", &ODE::OdeWrapper<2>::GetPropagatedTime)
        .def("SetStartTime", &ODE::OdeWrapper<2>::SetStartTime)
    ;

}


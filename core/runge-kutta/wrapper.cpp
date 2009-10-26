
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <rungekuttawrapper.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(librungekutta)
{
    scope* RungeKutta_RungeKuttaWrapper_1_scope = new scope(
    class_< RungeKutta::RungeKuttaWrapper<1>, boost::noncopyable >("RungeKutta_RungeKuttaWrapper_1", init<  >())
        .def_readwrite("ImTime", &RungeKutta::RungeKuttaWrapper<1>::ImTime)
        .def("ApplyConfigSection", &RungeKutta::RungeKuttaWrapper<1>::ApplyConfigSection)
        .def("Setup", &RungeKutta::RungeKuttaWrapper<1>::Setup)
        .def("AdvanceStep", &RungeKutta::RungeKuttaWrapper<1>::AdvanceStep)
    );

    enum_< RungeKutta::RungeKuttaWrapper<1>::IntegratorTypes >("IntegratorTypes")
        .value("IntegratorRK2", RungeKutta::RungeKuttaWrapper<1>::IntegratorRK2)
        .value("IntegratorRKCK", RungeKutta::RungeKuttaWrapper<1>::IntegratorRKCK)
        .value("IntegratorRK8PD", RungeKutta::RungeKuttaWrapper<1>::IntegratorRK8PD)
        .value("IntegratorRKF45", RungeKutta::RungeKuttaWrapper<1>::IntegratorRKF45)
        .value("IntegratorRK4", RungeKutta::RungeKuttaWrapper<1>::IntegratorRK4)
    ;

    delete RungeKutta_RungeKuttaWrapper_1_scope;

    scope* RungeKutta_RungeKuttaWrapper_2_scope = new scope(
    class_< RungeKutta::RungeKuttaWrapper<2>, boost::noncopyable >("RungeKutta_RungeKuttaWrapper_2", init<  >())
        .def_readwrite("ImTime", &RungeKutta::RungeKuttaWrapper<2>::ImTime)
        .def("ApplyConfigSection", &RungeKutta::RungeKuttaWrapper<2>::ApplyConfigSection)
        .def("Setup", &RungeKutta::RungeKuttaWrapper<2>::Setup)
        .def("AdvanceStep", &RungeKutta::RungeKuttaWrapper<2>::AdvanceStep)
    );

    enum_< RungeKutta::RungeKuttaWrapper<2>::IntegratorTypes >("IntegratorTypes")
        .value("IntegratorRK2", RungeKutta::RungeKuttaWrapper<2>::IntegratorRK2)
        .value("IntegratorRKCK", RungeKutta::RungeKuttaWrapper<2>::IntegratorRKCK)
        .value("IntegratorRK8PD", RungeKutta::RungeKuttaWrapper<2>::IntegratorRK8PD)
        .value("IntegratorRKF45", RungeKutta::RungeKuttaWrapper<2>::IntegratorRKF45)
        .value("IntegratorRK4", RungeKutta::RungeKuttaWrapper<2>::IntegratorRK4)
    ;

    delete RungeKutta_RungeKuttaWrapper_2_scope;

}


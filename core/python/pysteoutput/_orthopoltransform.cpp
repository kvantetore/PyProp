
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <transform/orthopol/orthopolpropagator.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_orthopoltransform()
{
    class_< OrthoPol::Propagator<1>, boost::noncopyable >("OrthoPolPropagator_1", no_init)
        .def("ApplyConfigSection", &OrthoPol::Propagator<1>::ApplyConfigSection)
        .def("Setup", &OrthoPol::Propagator<1>::Setup)
        .def("AdvanceStep", &OrthoPol::Propagator<1>::AdvanceStep)
        .def("ApplyDifferentiationMatrix", &OrthoPol::Propagator<1>::ApplyDifferentiationMatrix)
        .def("GetPropagationMatrix", &OrthoPol::Propagator<1>::GetPropagationMatrix)
        .def("GetDifferentiationMatrix", &OrthoPol::Propagator<1>::GetDifferentiationMatrix)
        .def("GetEigenvectors", &OrthoPol::Propagator<1>::GetEigenvectors)
        .def("GetEigenvalues", &OrthoPol::Propagator<1>::GetEigenvalues)
    ;

    class_< OrthoPol::Propagator<2>, boost::noncopyable >("OrthoPolPropagator_2", no_init)
        .def("ApplyConfigSection", &OrthoPol::Propagator<2>::ApplyConfigSection)
        .def("Setup", &OrthoPol::Propagator<2>::Setup)
        .def("AdvanceStep", &OrthoPol::Propagator<2>::AdvanceStep)
        .def("ApplyDifferentiationMatrix", &OrthoPol::Propagator<2>::ApplyDifferentiationMatrix)
        .def("GetPropagationMatrix", &OrthoPol::Propagator<2>::GetPropagationMatrix)
        .def("GetDifferentiationMatrix", &OrthoPol::Propagator<2>::GetDifferentiationMatrix)
        .def("GetEigenvectors", &OrthoPol::Propagator<2>::GetEigenvectors)
        .def("GetEigenvalues", &OrthoPol::Propagator<2>::GetEigenvalues)
    ;

    class_< OrthoPol::Propagator<3>, boost::noncopyable >("OrthoPolPropagator_3", no_init)
        .def("ApplyConfigSection", &OrthoPol::Propagator<3>::ApplyConfigSection)
        .def("Setup", &OrthoPol::Propagator<3>::Setup)
        .def("AdvanceStep", &OrthoPol::Propagator<3>::AdvanceStep)
        .def("ApplyDifferentiationMatrix", &OrthoPol::Propagator<3>::ApplyDifferentiationMatrix)
        .def("GetPropagationMatrix", &OrthoPol::Propagator<3>::GetPropagationMatrix)
        .def("GetDifferentiationMatrix", &OrthoPol::Propagator<3>::GetDifferentiationMatrix)
        .def("GetEigenvectors", &OrthoPol::Propagator<3>::GetEigenvectors)
        .def("GetEigenvalues", &OrthoPol::Propagator<3>::GetEigenvalues)
    ;

    class_< OrthoPol::Propagator<4>, boost::noncopyable >("OrthoPolPropagator_4", no_init)
        .def("ApplyConfigSection", &OrthoPol::Propagator<4>::ApplyConfigSection)
        .def("Setup", &OrthoPol::Propagator<4>::Setup)
        .def("AdvanceStep", &OrthoPol::Propagator<4>::AdvanceStep)
        .def("ApplyDifferentiationMatrix", &OrthoPol::Propagator<4>::ApplyDifferentiationMatrix)
        .def("GetPropagationMatrix", &OrthoPol::Propagator<4>::GetPropagationMatrix)
        .def("GetDifferentiationMatrix", &OrthoPol::Propagator<4>::GetDifferentiationMatrix)
        .def("GetEigenvectors", &OrthoPol::Propagator<4>::GetEigenvectors)
        .def("GetEigenvalues", &OrthoPol::Propagator<4>::GetEigenvalues)
    ;

}


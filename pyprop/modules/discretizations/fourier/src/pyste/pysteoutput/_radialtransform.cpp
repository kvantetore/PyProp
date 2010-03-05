
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <src/combinedpropagator/radialtransform.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_pyprop_modules_discretizations_fourier_src_pyste_radialtransform()
{
    class_< RadialTransform<1>, boost::noncopyable >("RadialTransform_1", init<  >())
        .def("TransformRank", &RadialTransform<1>::TransformRank)
        .def("ForwardTransform", &RadialTransform<1>::ForwardTransform)
        .def("InverseTransform", &RadialTransform<1>::InverseTransform)
    ;

    class_< RadialTransform<2>, boost::noncopyable >("RadialTransform_2", init<  >())
        .def("TransformRank", &RadialTransform<2>::TransformRank)
        .def("ForwardTransform", &RadialTransform<2>::ForwardTransform)
        .def("InverseTransform", &RadialTransform<2>::InverseTransform)
    ;

    class_< RadialTransform<3>, boost::noncopyable >("RadialTransform_3", init<  >())
        .def("TransformRank", &RadialTransform<3>::TransformRank)
        .def("ForwardTransform", &RadialTransform<3>::ForwardTransform)
        .def("InverseTransform", &RadialTransform<3>::InverseTransform)
    ;

}


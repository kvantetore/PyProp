
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/bsplinegridrepresentation.h>
#include <representation/bsplinerepresentation.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_core_modules_discretizations_bspline_pyste_bsplinerepresentation()
{
    scope* BSpline_BSplineRepresentation_scope = new scope(
    class_< BSpline::BSplineRepresentation, bases< Representation<1> >  >("BSplineRepresentation", init<  >())
        .def(init< const BSpline::BSplineRepresentation& >())
        .def("Copy", &BSpline::BSplineRepresentation::Copy)
        .def("GetFullShape", &BSpline::BSplineRepresentation::GetFullShape)
        .def("InnerProduct", &BSpline::BSplineRepresentation::InnerProduct)
        .def("GetGlobalWeights", &BSpline::BSplineRepresentation::GetGlobalWeights)
        .def("GetGlobalGrid", &BSpline::BSplineRepresentation::GetGlobalGrid)
        .def("GetGlobalOverlapMatrix", &BSpline::BSplineRepresentation::GetGlobalOverlapMatrix)
        .def("ApplyConfigSection", &BSpline::BSplineRepresentation::ApplyConfigSection)
        .def("SetupRepresentation", &BSpline::BSplineRepresentation::SetupRepresentation)
        .def("GetBSplineObject", &BSpline::BSplineRepresentation::GetBSplineObject)
    );
    register_ptr_to_python< boost::shared_ptr< BSpline::BSplineRepresentation > >();
    delete BSpline_BSplineRepresentation_scope;

    scope* BSpline_BSplineGridRepresentation_scope = new scope(
    class_< BSpline::BSplineGridRepresentation, bases< Representation<1> >  >("BSplineGridRepresentation", init<  >())
        .def(init< const BSpline::BSplineGridRepresentation& >())
        .def("Copy", &BSpline::BSplineGridRepresentation::Copy)
        .def("GetFullShape", &BSpline::BSplineGridRepresentation::GetFullShape)
        .def("InnerProduct", &BSpline::BSplineGridRepresentation::InnerProduct)
        .def("GetGlobalGrid", &BSpline::BSplineGridRepresentation::GetGlobalGrid)
        .def("GetGlobalWeights", &BSpline::BSplineGridRepresentation::GetGlobalWeights)
        .def("SetupRepresentation", &BSpline::BSplineGridRepresentation::SetupRepresentation)
        .def("GetBSplineObject", &BSpline::BSplineGridRepresentation::GetBSplineObject)
        .def("ApplyConfigSection", &BSpline::BSplineGridRepresentation::ApplyConfigSection)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< BSpline::BSplineGridRepresentation > >();
    delete BSpline_BSplineGridRepresentation_scope;

}


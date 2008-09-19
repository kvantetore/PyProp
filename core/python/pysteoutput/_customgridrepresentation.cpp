
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/customgridrepresentation.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_customgridrepresentation()
{
    scope* CustomGridRepresentation_scope = new scope(
    class_< CustomGridRepresentation, bases< Representation<1> >  >("CustomGridRepresentation", init<  >())
        .def(init< const CustomGridRepresentation& >())
        .def("Copy", &CustomGridRepresentation::Copy)
        .def("GetGlobalGrid", &CustomGridRepresentation::GetGlobalGrid)
        .def("GetGlobalWeights", &CustomGridRepresentation::GetGlobalWeights)
        .def("GetFullShape", &CustomGridRepresentation::GetFullShape)
        .def("ApplyConfigSection", &CustomGridRepresentation::ApplyConfigSection)
        .def("InnerProduct", &OrthogonalRepresentation::InnerProduct)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< CustomGridRepresentation > >();
    delete CustomGridRepresentation_scope;

}


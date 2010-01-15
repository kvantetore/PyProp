
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/orthopol/orthopolradialrepresentation.h>
#include <representation/orthopol/orthopolrange.h>
#include <transform/orthopol/orthopoltools.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_orthopolrepresentation()
{
    class_< OrthoPol::Parameter, boost::noncopyable >("OrthoPolParameter", no_init)
        .def_readwrite("Type", &OrthoPol::Parameter::Type)
        .def_readwrite("Scaling", &OrthoPol::Parameter::Scaling)
        .def_readwrite("HypersphericalRank", &OrthoPol::Parameter::HypersphericalRank)
    ;

    enum_< OrthoPol::PolynomialType >("OrthoPolType")
        .value("LaguerrePolynomial", OrthoPol::LaguerrePolynomial)
        .value("HermitePolynomial", OrthoPol::HermitePolynomial)
    ;

    class_< OrthoPol::OrthoPolRange, boost::noncopyable >("OrthoPolRange", init<  >())
        .def(init< const OrthoPol::Parameter&, int >())
        .def_readwrite("Count", &OrthoPol::OrthoPolRange::Count)
        .def_readwrite("Param", &OrthoPol::OrthoPolRange::Param)
        .def("GetGrid", &OrthoPol::OrthoPolRange::GetGrid, return_value_policy< return_by_value >())
        .def("GetWeights", &OrthoPol::OrthoPolRange::GetWeights, return_value_policy< return_by_value >())
        .def("Initialize", &OrthoPol::OrthoPolRange::Initialize)
    ;

    scope* OrthoPol_OrthoPolRadialRepresentation_scope = new scope(
    class_< OrthoPol::OrthoPolRadialRepresentation, bases< Representation<1> >  >("OrthoPolRadialRepresentation", init<  >())
        .def(init< const OrthoPol::OrthoPolRadialRepresentation& >())
        .def_readwrite("Range", &OrthoPol::OrthoPolRadialRepresentation::Range)
        .def("Copy", &OrthoPol::OrthoPolRadialRepresentation::Copy)
        .def("GetFullGrid", &OrthoPol::OrthoPolRadialRepresentation::GetFullGrid)
        .def("GetFullShape", &OrthoPol::OrthoPolRadialRepresentation::GetFullShape)
        .def("InnerProduct", &OrthoPol::OrthoPolRadialRepresentation::InnerProduct)
        .def("GetGlobalGrid", &OrthoPol::OrthoPolRadialRepresentation::GetGlobalGrid)
        .def("GetGlobalWeights", &OrthoPol::OrthoPolRadialRepresentation::GetGlobalWeights)
        .def("ApplyConfigSection", &OrthoPol::OrthoPolRadialRepresentation::ApplyConfigSection)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< OrthoPol::OrthoPolRadialRepresentation > >();
    delete OrthoPol_OrthoPolRadialRepresentation_scope;

}


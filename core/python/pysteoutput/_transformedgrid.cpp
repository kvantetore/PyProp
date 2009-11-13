
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/transformedgrid/transformedradialrepresentation.h>
#include <representation/transformedgrid/transformedrange.h>
#include <transform/transformedgrid/tools.h>
#include <transform/transformedgrid/transformedgridpropagator.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_transformedgrid()
{
    class_< TransformedGrid::Parameter, boost::noncopyable >("TransformedGridParameter", init<  >())
        .def_readwrite("Type", &TransformedGrid::Parameter::Type)
        .def_readwrite("Range", &TransformedGrid::Parameter::Range)
        .def_readwrite("Scaling", &TransformedGrid::Parameter::Scaling)
    ;

    enum_< TransformedGrid::TransformType >("TransformType")
        .value("TransformTypeLogarithmic", TransformedGrid::TransformTypeLogarithmic)
        .value("TransformTypeAlgebraic", TransformedGrid::TransformTypeAlgebraic)
        .value("TransformTypeTrigonometric", TransformedGrid::TransformTypeTrigonometric)
    ;

    enum_< TransformedGrid::TransformRange >("TransformRange")
        .value("TransformRangeRadial", TransformedGrid::TransformRangeRadial)
        .value("TransformRangeCartesian", TransformedGrid::TransformRangeCartesian)
    ;

    class_< TransformedRange, boost::noncopyable >("TransformedRange", init<  >())
        .def(init< const TransformedGrid::Parameter&, int >())
        .def_readwrite("Count", &TransformedRange::Count)
        .def_readwrite("Param", &TransformedRange::Param)
        .def("GetGrid", &TransformedRange::GetGrid, return_value_policy< return_by_value >())
        .def("GetWeights", &TransformedRange::GetWeights, return_value_policy< return_by_value >())
    ;

    class_< TransformedRadialRepresentation, bases< Representation<1> >  >("TransformedRadialRepresentation", init<  >())
        .def(init< const TransformedRadialRepresentation& >())
        .def_readwrite("Range", &TransformedRadialRepresentation::Range)
        .def("Copy", &TransformedRadialRepresentation::Copy)
        .def("GetFullShape", &TransformedRadialRepresentation::GetFullShape)
        .def("InnerProduct", &TransformedRadialRepresentation::InnerProduct)
        .def("GetGlobalGrid", &TransformedRadialRepresentation::GetGlobalGrid)
        .def("GetGlobalWeights", &TransformedRadialRepresentation::GetGlobalWeights)
        .def("ApplyConfigSection", &TransformedRadialRepresentation::ApplyConfigSection)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    ;

    class_< TransformedGrid::Propagator<1>, boost::noncopyable >("TransformedGridPropagator_1", no_init)
        .def("ApplyConfigSection", &TransformedGrid::Propagator<1>::ApplyConfigSection)
        .def("Setup", &TransformedGrid::Propagator<1>::Setup)
        .def("AdvanceStep", &TransformedGrid::Propagator<1>::AdvanceStep)
        .def("ApplyDifferentiationMatrix", &TransformedGrid::Propagator<1>::ApplyDifferentiationMatrix)
        .def("GetPropagationMatrix", &TransformedGrid::Propagator<1>::GetPropagationMatrix)
        .def("GetDifferentiationMatrix", &TransformedGrid::Propagator<1>::GetDifferentiationMatrix)
        .def("GetEigenvectors", &TransformedGrid::Propagator<1>::GetEigenvectors)
        .def("GetEigenvalues", &TransformedGrid::Propagator<1>::GetEigenvalues)
    ;

    class_< TransformedGrid::Propagator<2>, boost::noncopyable >("TransformedGridPropagator_2", no_init)
        .def("ApplyConfigSection", &TransformedGrid::Propagator<2>::ApplyConfigSection)
        .def("Setup", &TransformedGrid::Propagator<2>::Setup)
        .def("AdvanceStep", &TransformedGrid::Propagator<2>::AdvanceStep)
        .def("ApplyDifferentiationMatrix", &TransformedGrid::Propagator<2>::ApplyDifferentiationMatrix)
        .def("GetPropagationMatrix", &TransformedGrid::Propagator<2>::GetPropagationMatrix)
        .def("GetDifferentiationMatrix", &TransformedGrid::Propagator<2>::GetDifferentiationMatrix)
        .def("GetEigenvectors", &TransformedGrid::Propagator<2>::GetEigenvectors)
        .def("GetEigenvalues", &TransformedGrid::Propagator<2>::GetEigenvalues)
    ;

    class_< TransformedGrid::Propagator<3>, boost::noncopyable >("TransformedGridPropagator_3", no_init)
        .def("ApplyConfigSection", &TransformedGrid::Propagator<3>::ApplyConfigSection)
        .def("Setup", &TransformedGrid::Propagator<3>::Setup)
        .def("AdvanceStep", &TransformedGrid::Propagator<3>::AdvanceStep)
        .def("ApplyDifferentiationMatrix", &TransformedGrid::Propagator<3>::ApplyDifferentiationMatrix)
        .def("GetPropagationMatrix", &TransformedGrid::Propagator<3>::GetPropagationMatrix)
        .def("GetDifferentiationMatrix", &TransformedGrid::Propagator<3>::GetDifferentiationMatrix)
        .def("GetEigenvectors", &TransformedGrid::Propagator<3>::GetEigenvectors)
        .def("GetEigenvalues", &TransformedGrid::Propagator<3>::GetEigenvalues)
    ;

}


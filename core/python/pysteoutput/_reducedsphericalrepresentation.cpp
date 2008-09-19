
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/reducedspherical/lrange.h>
#include <representation/reducedspherical/reducedsphericalharmonicrepresentation.h>
#include <representation/reducedspherical/thetarange.h>
#include <representation/reducedspherical/thetarepresentation.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_reducedsphericalrepresentation()
{
    scope* ReducedSpherical_ReducedSphericalHarmonicRepresentation_scope = new scope(
    class_< ReducedSpherical::ReducedSphericalHarmonicRepresentation, bases< Representation<1> >  >("ReducedSphericalHarmonicRepresentation", init<  >())
        .def(init< const ReducedSpherical::ReducedSphericalHarmonicRepresentation& >())
        .def_readwrite("Range", &ReducedSpherical::ReducedSphericalHarmonicRepresentation::Range)
        .def("Copy", &ReducedSpherical::ReducedSphericalHarmonicRepresentation::Copy)
        .def("SetupRepresentation", &ReducedSpherical::ReducedSphericalHarmonicRepresentation::SetupRepresentation)
        .def("GetFullShape", &ReducedSpherical::ReducedSphericalHarmonicRepresentation::GetFullShape)
        .def("InnerProduct", &ReducedSpherical::ReducedSphericalHarmonicRepresentation::InnerProduct)
        .def("GetGlobalWeights", &ReducedSpherical::ReducedSphericalHarmonicRepresentation::GetGlobalWeights)
        .def("GetGlobalGrid", &ReducedSpherical::ReducedSphericalHarmonicRepresentation::GetGlobalGrid)
        .def("ApplyConfigSection", &ReducedSpherical::ReducedSphericalHarmonicRepresentation::ApplyConfigSection)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< ReducedSpherical::ReducedSphericalHarmonicRepresentation > >();
    delete ReducedSpherical_ReducedSphericalHarmonicRepresentation_scope;

    scope* ReducedSpherical_ThetaRepresentation_scope = new scope(
    class_< ReducedSpherical::ThetaRepresentation, bases< Representation<1> >  >("ThetaRepresentation", init<  >())
        .def(init< const ReducedSpherical::ThetaRepresentation& >())
        .def_readwrite("Range", &ReducedSpherical::ThetaRepresentation::Range)
        .def("Copy", &ReducedSpherical::ThetaRepresentation::Copy)
        .def("SetupRepresentation", &ReducedSpherical::ThetaRepresentation::SetupRepresentation)
        .def("GetFullShape", &ReducedSpherical::ThetaRepresentation::GetFullShape)
        .def("InnerProduct", &ReducedSpherical::ThetaRepresentation::InnerProduct)
        .def("GetGlobalGrid", &ReducedSpherical::ThetaRepresentation::GetGlobalGrid)
        .def("GetGlobalWeights", &ReducedSpherical::ThetaRepresentation::GetGlobalWeights)
        .def("ApplyConfigSection", &ReducedSpherical::ThetaRepresentation::ApplyConfigSection)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< ReducedSpherical::ThetaRepresentation > >();
    delete ReducedSpherical_ThetaRepresentation_scope;

    class_< ReducedSpherical::ThetaRange >("ThetaRange", init<  >())
        .def(init< const ReducedSpherical::ThetaRange& >())
        .def_readwrite("MaxL", &ReducedSpherical::ThetaRange::MaxL)
        .def("SetupRange", &ReducedSpherical::ThetaRange::SetupRange)
        .def("Count", &ReducedSpherical::ThetaRange::Count)
        .def("GetGrid", &ReducedSpherical::ThetaRange::GetGrid, return_value_policy< return_by_value >())
        .def("GetWeights", &ReducedSpherical::ThetaRange::GetWeights, return_value_policy< copy_const_reference >())
    ;

    class_< ReducedSpherical::LRange >("LRange", init<  >())
        .def(init< const ReducedSpherical::LRange& >())
        .def(init< int >())
        .def_readwrite("MaxL", &ReducedSpherical::LRange::MaxL)
        .def("Count", &ReducedSpherical::LRange::Count)
        .def("GetGrid", &ReducedSpherical::LRange::GetGrid, return_value_policy< copy_const_reference >())
        .def("GetWeights", &ReducedSpherical::LRange::GetWeights, return_value_policy< copy_const_reference >())
    ;

}


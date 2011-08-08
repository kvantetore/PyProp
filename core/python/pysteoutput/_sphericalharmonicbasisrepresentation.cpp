
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/sphericalbasis/lmbasisrange.h>
#include <representation/sphericalbasis/sphericalharmonicbasisrepresentation.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_sphericalharmonicbasisrepresentation()
{
    scope* SphericalBasis_SphericalHarmonicBasisRepresentation_scope = new scope(
    class_< SphericalBasis::SphericalHarmonicBasisRepresentation, bases< Representation<1> >  >("SphericalHarmonicBasisRepresentation", init<  >())
        .def(init< const SphericalBasis::SphericalHarmonicBasisRepresentation& >())
        .def_readwrite("Range", &SphericalBasis::SphericalHarmonicBasisRepresentation::Range)
        .def("Copy", &SphericalBasis::SphericalHarmonicBasisRepresentation::Copy)
        .def("GetFullShape", &SphericalBasis::SphericalHarmonicBasisRepresentation::GetFullShape)
        .def("GetGlobalWeights", &SphericalBasis::SphericalHarmonicBasisRepresentation::GetGlobalWeights)
        .def("GetGlobalGrid", &SphericalBasis::SphericalHarmonicBasisRepresentation::GetGlobalGrid)
        .def("ApplyConfigSection", &SphericalBasis::SphericalHarmonicBasisRepresentation::ApplyConfigSection)
        .def("InnerProduct", &OrthogonalRepresentation::InnerProduct)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< SphericalBasis::SphericalHarmonicBasisRepresentation > >();
    delete SphericalBasis_SphericalHarmonicBasisRepresentation_scope;

    scope* SphericalBasis_LmBasisRange_scope = new scope(
    class_< SphericalBasis::LmBasisRange >("LmBasisRange", init<  >())
        .def(init< const SphericalBasis::LmBasisRange& >())
        .def("Count", &SphericalBasis::LmBasisRange::Count)
        .def("GetGrid", &SphericalBasis::LmBasisRange::GetGrid, return_value_policy< copy_const_reference >())
        .def("GetWeights", &SphericalBasis::LmBasisRange::GetWeights, return_value_policy< copy_const_reference >())
        .def("GetGridIndex", &SphericalBasis::LmBasisRange::GetGridIndex)
        .def("IsGridIndex", &SphericalBasis::LmBasisRange::IsGridIndex)
        .def("GetLmIndex", &SphericalBasis::LmBasisRange::GetLmIndex)
        .def("BeginIndexList", &SphericalBasis::LmBasisRange::BeginIndexList)
        .def("AddIndex", &SphericalBasis::LmBasisRange::AddIndex)
        .def("EndIndexList", &SphericalBasis::LmBasisRange::EndIndexList)
    );
    register_ptr_to_python< boost::shared_ptr< SphericalBasis::LmBasisRange > >();
    delete SphericalBasis_LmBasisRange_scope;

    class_< SphericalBasis::LmIndex >("LmIndex", init<  >())
        .def(init< const SphericalBasis::LmIndex& >())
        .def(init< int, int >())
        .def_readwrite("l", &SphericalBasis::LmIndex::l)
        .def_readwrite("m", &SphericalBasis::LmIndex::m)
        .def("GetInt", &SphericalBasis::LmIndex::GetInt)
        .def( self == self )
    ;

    class_< SphericalBasis::ClebschGordan >("ClebschGordan", init<  >())
        .def(init< const SphericalBasis::ClebschGordan& >())
        .def("__call__", &SphericalBasis::ClebschGordan::operator ())
    ;

}


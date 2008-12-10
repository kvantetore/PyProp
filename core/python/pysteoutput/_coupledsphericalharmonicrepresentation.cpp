
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/coupledspherical/coupledrange.h>
#include <representation/coupledspherical/coupledsphericalharmonicrepresentation.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_coupledsphericalharmonicrepresentation()
{
    scope* CoupledSpherical_CoupledSphericalHarmonicRepresentation_scope = new scope(
    class_< CoupledSpherical::CoupledSphericalHarmonicRepresentation, bases< Representation<1> >  >("CoupledSphericalHarmonicRepresentation", init<  >())
        .def(init< const CoupledSpherical::CoupledSphericalHarmonicRepresentation& >())
        .def_readwrite("Range", &CoupledSpherical::CoupledSphericalHarmonicRepresentation::Range)
        .def("Copy", &CoupledSpherical::CoupledSphericalHarmonicRepresentation::Copy)
        .def("GetFullShape", &CoupledSpherical::CoupledSphericalHarmonicRepresentation::GetFullShape)
        .def("GetGlobalWeights", &CoupledSpherical::CoupledSphericalHarmonicRepresentation::GetGlobalWeights)
        .def("GetGlobalGrid", &CoupledSpherical::CoupledSphericalHarmonicRepresentation::GetGlobalGrid)
        .def("ApplyConfigSection", &CoupledSpherical::CoupledSphericalHarmonicRepresentation::ApplyConfigSection)
        .def("InnerProduct", &OrthogonalRepresentation::InnerProduct)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< CoupledSpherical::CoupledSphericalHarmonicRepresentation > >();
    delete CoupledSpherical_CoupledSphericalHarmonicRepresentation_scope;

    scope* CoupledSpherical_CoupledRange_scope = new scope(
    class_< CoupledSpherical::CoupledRange, boost::noncopyable >("CoupledRange", init<  >())
        .def("Count", &CoupledSpherical::CoupledRange::Count)
        .def("GetGrid", &CoupledSpherical::CoupledRange::GetGrid, return_value_policy< copy_const_reference >())
        .def("GetWeights", &CoupledSpherical::CoupledRange::GetWeights, return_value_policy< copy_const_reference >())
        .def("GetGridIndex", &CoupledSpherical::CoupledRange::GetGridIndex)
        .def("IsGridIndex", &CoupledSpherical::CoupledRange::IsGridIndex)
        .def("GetCoupledIndex", &CoupledSpherical::CoupledRange::GetCoupledIndex)
        .def("BeginIndexList", &CoupledSpherical::CoupledRange::BeginIndexList)
        .def("AddIndex", &CoupledSpherical::CoupledRange::AddIndex)
        .def("EndIndexList", &CoupledSpherical::CoupledRange::EndIndexList)
    );
    register_ptr_to_python< boost::shared_ptr< CoupledSpherical::CoupledRange > >();
    delete CoupledSpherical_CoupledRange_scope;

    class_< CoupledSpherical::CoupledIndex >("CoupledIndex", init<  >())
        .def(init< const CoupledSpherical::CoupledIndex& >())
        .def(init< int, int, int, int >())
        .def_readwrite("l1", &CoupledSpherical::CoupledIndex::l1)
        .def_readwrite("l2", &CoupledSpherical::CoupledIndex::l2)
        .def_readwrite("L", &CoupledSpherical::CoupledIndex::L)
        .def_readwrite("M", &CoupledSpherical::CoupledIndex::M)
        .def("GetInt", &CoupledSpherical::CoupledIndex::GetInt)
        .def( self == self )
    ;

}


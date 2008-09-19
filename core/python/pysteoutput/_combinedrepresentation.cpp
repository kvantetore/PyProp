
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/combinedrepresentation.h>
#include <representation/spherical/angularrepresentation.h>
#include <representation/spherical/sphericalharmonicrepresentation.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

void (CombinedRepresentation<1>::*CombinedRepresentation_1___MultiplyOverlapWavefunction_1___Wavefunction_1___int)(Wavefunction<1>&, Wavefunction<1>&, int)  = &CombinedRepresentation<1>::MultiplyOverlap;

void (CombinedRepresentation<1>::*CombinedRepresentation_1___MultiplyOverlapWavefunction_1__)(Wavefunction<1>&)  = &CombinedRepresentation<1>::MultiplyOverlap;

void (CombinedRepresentation<2>::*CombinedRepresentation_2___MultiplyOverlapWavefunction_2___Wavefunction_2___int)(Wavefunction<2>&, Wavefunction<2>&, int)  = &CombinedRepresentation<2>::MultiplyOverlap;

void (CombinedRepresentation<2>::*CombinedRepresentation_2___MultiplyOverlapWavefunction_2__)(Wavefunction<2>&)  = &CombinedRepresentation<2>::MultiplyOverlap;

void (CombinedRepresentation<3>::*CombinedRepresentation_3___MultiplyOverlapWavefunction_3___Wavefunction_3___int)(Wavefunction<3>&, Wavefunction<3>&, int)  = &CombinedRepresentation<3>::MultiplyOverlap;

void (CombinedRepresentation<3>::*CombinedRepresentation_3___MultiplyOverlapWavefunction_3__)(Wavefunction<3>&)  = &CombinedRepresentation<3>::MultiplyOverlap;

void (CombinedRepresentation<4>::*CombinedRepresentation_4___MultiplyOverlapWavefunction_4___Wavefunction_4___int)(Wavefunction<4>&, Wavefunction<4>&, int)  = &CombinedRepresentation<4>::MultiplyOverlap;

void (CombinedRepresentation<4>::*CombinedRepresentation_4___MultiplyOverlapWavefunction_4__)(Wavefunction<4>&)  = &CombinedRepresentation<4>::MultiplyOverlap;


}// namespace 


// Module ======================================================================
void Export_python_combinedrepresentation()
{
    class_< CombinedRepresentation<1>, bases< Representation<1> >  >("CombinedRepresentation_1", init<  >())
        .def(init< const CombinedRepresentation<1>& >())
        .def_readwrite("Algorithm", &CombinedRepresentation<1>::Algorithm)
        .def("Copy", &CombinedRepresentation<1>::Copy)
        .def("GetRepresentation", &CombinedRepresentation<1>::GetRepresentation)
        .def("SetRepresentation", &CombinedRepresentation<1>::SetRepresentation)
        .def("GetFullShape", &CombinedRepresentation<1>::GetFullShape)
        .def("InnerProduct", &CombinedRepresentation<1>::InnerProduct)
        .def("GetGlobalGrid", &CombinedRepresentation<1>::GetGlobalGrid)
        .def("GetGlobalWeights", &CombinedRepresentation<1>::GetGlobalWeights)
        .def("ApplyConfigSection", &CombinedRepresentation<1>::ApplyConfigSection)
        .def("GetGlobalOverlapMatrix", &CombinedRepresentation<1>::GetGlobalOverlapMatrix)
        .def("MultiplyOverlap", CombinedRepresentation_1___MultiplyOverlapWavefunction_1___Wavefunction_1___int)
        .def("MultiplyOverlap", CombinedRepresentation_1___MultiplyOverlapWavefunction_1__)
        .def("SolveOverlap", &CombinedRepresentation<1>::SolveOverlap)
        .def("MultiplySqrtOverlap", &CombinedRepresentation<1>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &CombinedRepresentation<1>::SolveSqrtOverlap)
    ;

    class_< CombinedRepresentation<2>, bases< Representation<2> >  >("CombinedRepresentation_2", init<  >())
        .def(init< const CombinedRepresentation<2>& >())
        .def_readwrite("Algorithm", &CombinedRepresentation<2>::Algorithm)
        .def("Copy", &CombinedRepresentation<2>::Copy)
        .def("GetRepresentation", &CombinedRepresentation<2>::GetRepresentation)
        .def("SetRepresentation", &CombinedRepresentation<2>::SetRepresentation)
        .def("GetFullShape", &CombinedRepresentation<2>::GetFullShape)
        .def("InnerProduct", &CombinedRepresentation<2>::InnerProduct)
        .def("GetGlobalGrid", &CombinedRepresentation<2>::GetGlobalGrid)
        .def("GetGlobalWeights", &CombinedRepresentation<2>::GetGlobalWeights)
        .def("ApplyConfigSection", &CombinedRepresentation<2>::ApplyConfigSection)
        .def("GetGlobalOverlapMatrix", &CombinedRepresentation<2>::GetGlobalOverlapMatrix)
        .def("MultiplyOverlap", CombinedRepresentation_2___MultiplyOverlapWavefunction_2___Wavefunction_2___int)
        .def("MultiplyOverlap", CombinedRepresentation_2___MultiplyOverlapWavefunction_2__)
        .def("SolveOverlap", &CombinedRepresentation<2>::SolveOverlap)
        .def("MultiplySqrtOverlap", &CombinedRepresentation<2>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &CombinedRepresentation<2>::SolveSqrtOverlap)
    ;

    class_< CombinedRepresentation<3>, bases< Representation<3> >  >("CombinedRepresentation_3", init<  >())
        .def(init< const CombinedRepresentation<3>& >())
        .def_readwrite("Algorithm", &CombinedRepresentation<3>::Algorithm)
        .def("Copy", &CombinedRepresentation<3>::Copy)
        .def("GetRepresentation", &CombinedRepresentation<3>::GetRepresentation)
        .def("SetRepresentation", &CombinedRepresentation<3>::SetRepresentation)
        .def("GetFullShape", &CombinedRepresentation<3>::GetFullShape)
        .def("InnerProduct", &CombinedRepresentation<3>::InnerProduct)
        .def("GetGlobalGrid", &CombinedRepresentation<3>::GetGlobalGrid)
        .def("GetGlobalWeights", &CombinedRepresentation<3>::GetGlobalWeights)
        .def("ApplyConfigSection", &CombinedRepresentation<3>::ApplyConfigSection)
        .def("GetGlobalOverlapMatrix", &CombinedRepresentation<3>::GetGlobalOverlapMatrix)
        .def("MultiplyOverlap", CombinedRepresentation_3___MultiplyOverlapWavefunction_3___Wavefunction_3___int)
        .def("MultiplyOverlap", CombinedRepresentation_3___MultiplyOverlapWavefunction_3__)
        .def("SolveOverlap", &CombinedRepresentation<3>::SolveOverlap)
        .def("MultiplySqrtOverlap", &CombinedRepresentation<3>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &CombinedRepresentation<3>::SolveSqrtOverlap)
    ;

    class_< CombinedRepresentation<4>, bases< Representation<4> >  >("CombinedRepresentation_4", init<  >())
        .def(init< const CombinedRepresentation<4>& >())
        .def_readwrite("Algorithm", &CombinedRepresentation<4>::Algorithm)
        .def("Copy", &CombinedRepresentation<4>::Copy)
        .def("GetRepresentation", &CombinedRepresentation<4>::GetRepresentation)
        .def("SetRepresentation", &CombinedRepresentation<4>::SetRepresentation)
        .def("GetFullShape", &CombinedRepresentation<4>::GetFullShape)
        .def("InnerProduct", &CombinedRepresentation<4>::InnerProduct)
        .def("GetGlobalGrid", &CombinedRepresentation<4>::GetGlobalGrid)
        .def("GetGlobalWeights", &CombinedRepresentation<4>::GetGlobalWeights)
        .def("ApplyConfigSection", &CombinedRepresentation<4>::ApplyConfigSection)
        .def("GetGlobalOverlapMatrix", &CombinedRepresentation<4>::GetGlobalOverlapMatrix)
        .def("MultiplyOverlap", CombinedRepresentation_4___MultiplyOverlapWavefunction_4___Wavefunction_4___int)
        .def("MultiplyOverlap", CombinedRepresentation_4___MultiplyOverlapWavefunction_4__)
        .def("SolveOverlap", &CombinedRepresentation<4>::SolveOverlap)
        .def("MultiplySqrtOverlap", &CombinedRepresentation<4>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &CombinedRepresentation<4>::SolveSqrtOverlap)
    ;

    scope* AngularRepresentation_scope = new scope(
    class_< AngularRepresentation, bases< Representation<1> >  >("AngularRepresentation", init<  >())
        .def(init< const AngularRepresentation& >())
        .def_readwrite("Range", &AngularRepresentation::Range)
        .def("Copy", &AngularRepresentation::Copy)
        .def("SetupRepresentation", &AngularRepresentation::SetupRepresentation)
        .def("GetFullShape", &AngularRepresentation::GetFullShape)
        .def("InnerProduct", &AngularRepresentation::InnerProduct)
        .def("GetGlobalGrid", &AngularRepresentation::GetGlobalGrid)
        .def("GetGlobalWeights", &AngularRepresentation::GetGlobalWeights)
        .def("GetGlobalExpandedGrid", &AngularRepresentation::GetGlobalExpandedGrid)
        .def("ApplyConfigSection", &AngularRepresentation::ApplyConfigSection)
        .def("GetLocalExpandedGrid", &CompressedRepresentation::GetLocalExpandedGrid)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< AngularRepresentation > >();
    delete AngularRepresentation_scope;

    scope* SphericalHarmonicRepresentation_scope = new scope(
    class_< SphericalHarmonicRepresentation, bases< Representation<1> >  >("SphericalHarmonicRepresentation", init<  >())
        .def(init< const SphericalHarmonicRepresentation& >())
        .def_readwrite("Range", &SphericalHarmonicRepresentation::Range)
        .def("Copy", &SphericalHarmonicRepresentation::Copy)
        .def("SetupRepresentation", &SphericalHarmonicRepresentation::SetupRepresentation)
        .def("GetFullShape", &SphericalHarmonicRepresentation::GetFullShape)
        .def("InnerProduct", &SphericalHarmonicRepresentation::InnerProduct)
        .def("GetGlobalWeights", &SphericalHarmonicRepresentation::GetGlobalWeights)
        .def("GetGlobalGrid", &SphericalHarmonicRepresentation::GetGlobalGrid)
        .def("GetGlobalExpandedGrid", &SphericalHarmonicRepresentation::GetGlobalExpandedGrid)
        .def("ApplyConfigSection", &SphericalHarmonicRepresentation::ApplyConfigSection)
        .def("GetLocalExpandedGrid", &CompressedRepresentation::GetLocalExpandedGrid)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< SphericalHarmonicRepresentation > >();
    delete SphericalHarmonicRepresentation_scope;

}


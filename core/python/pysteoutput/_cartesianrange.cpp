
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/cartesianrange.h>
#include <representation/cartesianrepresentation.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

void (CartesianRepresentation<1>::*CartesianRepresentation_1___MultiplyOverlapWavefunction_1___Wavefunction_1___int)(Wavefunction<1>&, Wavefunction<1>&, int)  = &CartesianRepresentation<1>::MultiplyOverlap;

void (CartesianRepresentation<1>::*CartesianRepresentation_1___MultiplyOverlapWavefunction_1__)(Wavefunction<1>&)  = &CartesianRepresentation<1>::MultiplyOverlap;

void (CartesianRepresentation<2>::*CartesianRepresentation_2___MultiplyOverlapWavefunction_2___Wavefunction_2___int)(Wavefunction<2>&, Wavefunction<2>&, int)  = &CartesianRepresentation<2>::MultiplyOverlap;

void (CartesianRepresentation<2>::*CartesianRepresentation_2___MultiplyOverlapWavefunction_2__)(Wavefunction<2>&)  = &CartesianRepresentation<2>::MultiplyOverlap;

void (CartesianRepresentation<3>::*CartesianRepresentation_3___MultiplyOverlapWavefunction_3___Wavefunction_3___int)(Wavefunction<3>&, Wavefunction<3>&, int)  = &CartesianRepresentation<3>::MultiplyOverlap;

void (CartesianRepresentation<3>::*CartesianRepresentation_3___MultiplyOverlapWavefunction_3__)(Wavefunction<3>&)  = &CartesianRepresentation<3>::MultiplyOverlap;

void (CartesianRepresentation<4>::*CartesianRepresentation_4___MultiplyOverlapWavefunction_4___Wavefunction_4___int)(Wavefunction<4>&, Wavefunction<4>&, int)  = &CartesianRepresentation<4>::MultiplyOverlap;

void (CartesianRepresentation<4>::*CartesianRepresentation_4___MultiplyOverlapWavefunction_4__)(Wavefunction<4>&)  = &CartesianRepresentation<4>::MultiplyOverlap;


}// namespace 


// Module ======================================================================
void Export_python_cartesianrange()
{
    scope* CartesianRepresentation_1_scope = new scope(
    class_< CartesianRepresentation<1>, bases< Representation<1> >  >("CartesianRepresentation_1", init<  >())
        .def(init< const CartesianRepresentation<1>& >())
        .def(init< CartesianRange& >())
        .def(init< blitz::TinyVector<CartesianRange,1>& >())
        .def_readwrite("Range", &CartesianRepresentation<1>::Range)
        .def("Copy", &CartesianRepresentation<1>::Copy)
        .def("GetRange", &CartesianRepresentation<1>::GetRange, return_value_policy< copy_const_reference >())
        .def("GetGlobalGrid", &CartesianRepresentation<1>::GetGlobalGrid)
        .def("GetGlobalWeights", &CartesianRepresentation<1>::GetGlobalWeights)
        .def("GetScalarWeight", &CartesianRepresentation<1>::GetScalarWeight)
        .def("GetFullShape", &CartesianRepresentation<1>::GetFullShape)
        .def("InnerProduct", &CartesianRepresentation<1>::InnerProduct)
        .def("ApplyConfigSection", &CartesianRepresentation<1>::ApplyConfigSection)
        .def("MultiplyOverlap", CartesianRepresentation_1___MultiplyOverlapWavefunction_1___Wavefunction_1___int)
        .def("MultiplyOverlap", CartesianRepresentation_1___MultiplyOverlapWavefunction_1__)
        .def("SolveOverlap", &CartesianRepresentation<1>::SolveOverlap)
        .def("MultiplySqrtOverlap", &CartesianRepresentation<1>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &CartesianRepresentation<1>::SolveSqrtOverlap)
        .def("GetGlobalOverlapMatrix", &CartesianRepresentation<1>::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< CartesianRepresentation<1> > >();
    delete CartesianRepresentation_1_scope;

    scope* CartesianRepresentation_2_scope = new scope(
    class_< CartesianRepresentation<2>, bases< Representation<2> >  >("CartesianRepresentation_2", init<  >())
        .def(init< const CartesianRepresentation<2>& >())
        .def(init< CartesianRange& >())
        .def(init< blitz::TinyVector<CartesianRange,2>& >())
        .def_readwrite("Range", &CartesianRepresentation<2>::Range)
        .def("Copy", &CartesianRepresentation<2>::Copy)
        .def("GetRange", &CartesianRepresentation<2>::GetRange, return_value_policy< copy_const_reference >())
        .def("GetGlobalGrid", &CartesianRepresentation<2>::GetGlobalGrid)
        .def("GetGlobalWeights", &CartesianRepresentation<2>::GetGlobalWeights)
        .def("GetScalarWeight", &CartesianRepresentation<2>::GetScalarWeight)
        .def("GetFullShape", &CartesianRepresentation<2>::GetFullShape)
        .def("InnerProduct", &CartesianRepresentation<2>::InnerProduct)
        .def("ApplyConfigSection", &CartesianRepresentation<2>::ApplyConfigSection)
        .def("MultiplyOverlap", CartesianRepresentation_2___MultiplyOverlapWavefunction_2___Wavefunction_2___int)
        .def("MultiplyOverlap", CartesianRepresentation_2___MultiplyOverlapWavefunction_2__)
        .def("SolveOverlap", &CartesianRepresentation<2>::SolveOverlap)
        .def("MultiplySqrtOverlap", &CartesianRepresentation<2>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &CartesianRepresentation<2>::SolveSqrtOverlap)
        .def("GetGlobalOverlapMatrix", &CartesianRepresentation<2>::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< CartesianRepresentation<2> > >();
    delete CartesianRepresentation_2_scope;

    scope* CartesianRepresentation_3_scope = new scope(
    class_< CartesianRepresentation<3>, bases< Representation<3> >  >("CartesianRepresentation_3", init<  >())
        .def(init< const CartesianRepresentation<3>& >())
        .def(init< CartesianRange& >())
        .def(init< blitz::TinyVector<CartesianRange,3>& >())
        .def_readwrite("Range", &CartesianRepresentation<3>::Range)
        .def("Copy", &CartesianRepresentation<3>::Copy)
        .def("GetRange", &CartesianRepresentation<3>::GetRange, return_value_policy< copy_const_reference >())
        .def("GetGlobalGrid", &CartesianRepresentation<3>::GetGlobalGrid)
        .def("GetGlobalWeights", &CartesianRepresentation<3>::GetGlobalWeights)
        .def("GetScalarWeight", &CartesianRepresentation<3>::GetScalarWeight)
        .def("GetFullShape", &CartesianRepresentation<3>::GetFullShape)
        .def("InnerProduct", &CartesianRepresentation<3>::InnerProduct)
        .def("ApplyConfigSection", &CartesianRepresentation<3>::ApplyConfigSection)
        .def("MultiplyOverlap", CartesianRepresentation_3___MultiplyOverlapWavefunction_3___Wavefunction_3___int)
        .def("MultiplyOverlap", CartesianRepresentation_3___MultiplyOverlapWavefunction_3__)
        .def("SolveOverlap", &CartesianRepresentation<3>::SolveOverlap)
        .def("MultiplySqrtOverlap", &CartesianRepresentation<3>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &CartesianRepresentation<3>::SolveSqrtOverlap)
        .def("GetGlobalOverlapMatrix", &CartesianRepresentation<3>::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< CartesianRepresentation<3> > >();
    delete CartesianRepresentation_3_scope;

    scope* CartesianRepresentation_4_scope = new scope(
    class_< CartesianRepresentation<4>, bases< Representation<4> >  >("CartesianRepresentation_4", init<  >())
        .def(init< const CartesianRepresentation<4>& >())
        .def(init< CartesianRange& >())
        .def(init< blitz::TinyVector<CartesianRange,4>& >())
        .def_readwrite("Range", &CartesianRepresentation<4>::Range)
        .def("Copy", &CartesianRepresentation<4>::Copy)
        .def("GetRange", &CartesianRepresentation<4>::GetRange, return_value_policy< copy_const_reference >())
        .def("GetGlobalGrid", &CartesianRepresentation<4>::GetGlobalGrid)
        .def("GetGlobalWeights", &CartesianRepresentation<4>::GetGlobalWeights)
        .def("GetScalarWeight", &CartesianRepresentation<4>::GetScalarWeight)
        .def("GetFullShape", &CartesianRepresentation<4>::GetFullShape)
        .def("InnerProduct", &CartesianRepresentation<4>::InnerProduct)
        .def("ApplyConfigSection", &CartesianRepresentation<4>::ApplyConfigSection)
        .def("MultiplyOverlap", CartesianRepresentation_4___MultiplyOverlapWavefunction_4___Wavefunction_4___int)
        .def("MultiplyOverlap", CartesianRepresentation_4___MultiplyOverlapWavefunction_4__)
        .def("SolveOverlap", &CartesianRepresentation<4>::SolveOverlap)
        .def("MultiplySqrtOverlap", &CartesianRepresentation<4>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &CartesianRepresentation<4>::SolveSqrtOverlap)
        .def("GetGlobalOverlapMatrix", &CartesianRepresentation<4>::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< CartesianRepresentation<4> > >();
    delete CartesianRepresentation_4_scope;

    class_< CartesianRange >("CartesianRange", init<  >())
        .def(init< double, double, int, optional< bool > >())
        .def(init< const CartesianRange& >())
        .def_readwrite("Min", &CartesianRange::Min)
        .def_readwrite("Max", &CartesianRange::Max)
        .def_readwrite("Dx", &CartesianRange::Dx)
        .def_readwrite("Count", &CartesianRange::Count)
        .def_readwrite("TranslatedGrid", &CartesianRange::TranslatedGrid)
        .def("GetGrid", &CartesianRange::GetGrid, return_value_policy< return_by_value >())
        .def("GetWeights", &CartesianRange::GetWeights, return_value_policy< return_by_value >())
        .def("GetPosition", &CartesianRange::GetPosition)
        .def("GetIndex", &CartesianRange::GetIndex)
        .def("GetOverlapMatrix", &CartesianRange::GetOverlapMatrix)
    ;

}


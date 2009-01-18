
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/representation.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

void (Representation<1>::*Representation_1___MultiplyIntegrationWeightsWavefunction_1___Wavefunction_1___int)(Wavefunction<1>&, Wavefunction<1>&, int)  = &Representation<1>::MultiplyIntegrationWeights;

void (Representation<1>::*Representation_1___MultiplyIntegrationWeightsWavefunction_1___Wavefunction_1__)(Wavefunction<1>&, Wavefunction<1>&)  = &Representation<1>::MultiplyIntegrationWeights;

void (Representation<1>::*Representation_1___MultiplyIntegrationWeightsWavefunction_1___int)(Wavefunction<1>&, int)  = &Representation<1>::MultiplyIntegrationWeights;

void (Representation<1>::*Representation_1___MultiplyIntegrationWeightsWavefunction_1__)(Wavefunction<1>&)  = &Representation<1>::MultiplyIntegrationWeights;

void (Representation<1>::*Representation_1___MultiplyOverlapWavefunction_1___Wavefunction_1___int)(Wavefunction<1>&, Wavefunction<1>&, int)  = &Representation<1>::MultiplyOverlap;

void (Representation<1>::*Representation_1___MultiplyOverlapstd__complex_double__Wavefunction_1___std__complex_double__Wavefunction_1___int)(std::complex<double>, Wavefunction<1>&, std::complex<double>, Wavefunction<1>&, int)  = &Representation<1>::MultiplyOverlap;

void (Representation<1>::*Representation_1___MultiplyOverlapWavefunction_1__)(Wavefunction<1>&)  = &Representation<1>::MultiplyOverlap;

void (Representation<2>::*Representation_2___MultiplyIntegrationWeightsWavefunction_2___Wavefunction_2___int)(Wavefunction<2>&, Wavefunction<2>&, int)  = &Representation<2>::MultiplyIntegrationWeights;

void (Representation<2>::*Representation_2___MultiplyIntegrationWeightsWavefunction_2___Wavefunction_2__)(Wavefunction<2>&, Wavefunction<2>&)  = &Representation<2>::MultiplyIntegrationWeights;

void (Representation<2>::*Representation_2___MultiplyIntegrationWeightsWavefunction_2___int)(Wavefunction<2>&, int)  = &Representation<2>::MultiplyIntegrationWeights;

void (Representation<2>::*Representation_2___MultiplyIntegrationWeightsWavefunction_2__)(Wavefunction<2>&)  = &Representation<2>::MultiplyIntegrationWeights;

void (Representation<2>::*Representation_2___MultiplyOverlapWavefunction_2___Wavefunction_2___int)(Wavefunction<2>&, Wavefunction<2>&, int)  = &Representation<2>::MultiplyOverlap;

void (Representation<2>::*Representation_2___MultiplyOverlapstd__complex_double__Wavefunction_2___std__complex_double__Wavefunction_2___int)(std::complex<double>, Wavefunction<2>&, std::complex<double>, Wavefunction<2>&, int)  = &Representation<2>::MultiplyOverlap;

void (Representation<2>::*Representation_2___MultiplyOverlapWavefunction_2__)(Wavefunction<2>&)  = &Representation<2>::MultiplyOverlap;

void (Representation<3>::*Representation_3___MultiplyIntegrationWeightsWavefunction_3___Wavefunction_3___int)(Wavefunction<3>&, Wavefunction<3>&, int)  = &Representation<3>::MultiplyIntegrationWeights;

void (Representation<3>::*Representation_3___MultiplyIntegrationWeightsWavefunction_3___Wavefunction_3__)(Wavefunction<3>&, Wavefunction<3>&)  = &Representation<3>::MultiplyIntegrationWeights;

void (Representation<3>::*Representation_3___MultiplyIntegrationWeightsWavefunction_3___int)(Wavefunction<3>&, int)  = &Representation<3>::MultiplyIntegrationWeights;

void (Representation<3>::*Representation_3___MultiplyIntegrationWeightsWavefunction_3__)(Wavefunction<3>&)  = &Representation<3>::MultiplyIntegrationWeights;

void (Representation<3>::*Representation_3___MultiplyOverlapWavefunction_3___Wavefunction_3___int)(Wavefunction<3>&, Wavefunction<3>&, int)  = &Representation<3>::MultiplyOverlap;

void (Representation<3>::*Representation_3___MultiplyOverlapstd__complex_double__Wavefunction_3___std__complex_double__Wavefunction_3___int)(std::complex<double>, Wavefunction<3>&, std::complex<double>, Wavefunction<3>&, int)  = &Representation<3>::MultiplyOverlap;

void (Representation<3>::*Representation_3___MultiplyOverlapWavefunction_3__)(Wavefunction<3>&)  = &Representation<3>::MultiplyOverlap;

void (Representation<4>::*Representation_4___MultiplyIntegrationWeightsWavefunction_4___Wavefunction_4___int)(Wavefunction<4>&, Wavefunction<4>&, int)  = &Representation<4>::MultiplyIntegrationWeights;

void (Representation<4>::*Representation_4___MultiplyIntegrationWeightsWavefunction_4___Wavefunction_4__)(Wavefunction<4>&, Wavefunction<4>&)  = &Representation<4>::MultiplyIntegrationWeights;

void (Representation<4>::*Representation_4___MultiplyIntegrationWeightsWavefunction_4___int)(Wavefunction<4>&, int)  = &Representation<4>::MultiplyIntegrationWeights;

void (Representation<4>::*Representation_4___MultiplyIntegrationWeightsWavefunction_4__)(Wavefunction<4>&)  = &Representation<4>::MultiplyIntegrationWeights;

void (Representation<4>::*Representation_4___MultiplyOverlapWavefunction_4___Wavefunction_4___int)(Wavefunction<4>&, Wavefunction<4>&, int)  = &Representation<4>::MultiplyOverlap;

void (Representation<4>::*Representation_4___MultiplyOverlapstd__complex_double__Wavefunction_4___std__complex_double__Wavefunction_4___int)(std::complex<double>, Wavefunction<4>&, std::complex<double>, Wavefunction<4>&, int)  = &Representation<4>::MultiplyOverlap;

void (Representation<4>::*Representation_4___MultiplyOverlapWavefunction_4__)(Wavefunction<4>&)  = &Representation<4>::MultiplyOverlap;


}// namespace 


// Module ======================================================================
void Export_python_representation()
{
    scope* Representation_1_scope = new scope(
    class_< Representation<1>, boost::noncopyable >("Representation_1", no_init)
        .def("SetDistributedModel", &Representation<1>::SetDistributedModel)
        .def("GetDistributedModel", &Representation<1>::GetDistributedModel)
        .def("GetInitialShape", &Representation<1>::GetInitialShape)
        .def("GetBaseRank", &Representation<1>::GetBaseRank)
        .def("SetBaseRank", &Representation<1>::SetBaseRank)
        .def("GetId", &Representation<1>::GetId)
        .def("GetLocalGrid", &Representation<1>::GetLocalGrid)
        .def("GetLocalWeights", &Representation<1>::GetLocalWeights)
        .def("MultiplyIntegrationWeights", Representation_1___MultiplyIntegrationWeightsWavefunction_1___Wavefunction_1___int)
        .def("MultiplyIntegrationWeights", Representation_1___MultiplyIntegrationWeightsWavefunction_1___Wavefunction_1__)
        .def("MultiplyIntegrationWeights", Representation_1___MultiplyIntegrationWeightsWavefunction_1___int)
        .def("MultiplyIntegrationWeights", Representation_1___MultiplyIntegrationWeightsWavefunction_1__)
        .def("GetGlobalOverlapMatrix", &Representation<1>::GetGlobalOverlapMatrix)
        .def("MultiplyOverlap", Representation_1___MultiplyOverlapWavefunction_1___Wavefunction_1___int)
        .def("MultiplyOverlap", Representation_1___MultiplyOverlapstd__complex_double__Wavefunction_1___std__complex_double__Wavefunction_1___int)
        .def("MultiplyOverlap", Representation_1___MultiplyOverlapWavefunction_1__)
        .def("SolveOverlap", &Representation<1>::SolveOverlap)
        .def("MultiplySqrtOverlap", &Representation<1>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &Representation<1>::SolveSqrtOverlap)
        .def("IsOrthogonalBasis", &Representation<1>::IsOrthogonalBasis)
        .def("GetFullShape", &Representation<1>::GetFullShape)
        .def("InnerProduct", &Representation<1>::InnerProduct)
        .def("GetGlobalWeights", &Representation<1>::GetGlobalWeights)
        .def("GetGlobalGrid", &Representation<1>::GetGlobalGrid)
        .def("ApplyConfigSection", &Representation<1>::ApplyConfigSection)
        .def("Copy", &Representation<1>::Copy)
    );
    register_ptr_to_python< boost::shared_ptr< Representation<1> > >();
    delete Representation_1_scope;

    scope* Representation_2_scope = new scope(
    class_< Representation<2>, boost::noncopyable >("Representation_2", no_init)
        .def("SetDistributedModel", &Representation<2>::SetDistributedModel)
        .def("GetDistributedModel", &Representation<2>::GetDistributedModel)
        .def("GetInitialShape", &Representation<2>::GetInitialShape)
        .def("GetBaseRank", &Representation<2>::GetBaseRank)
        .def("SetBaseRank", &Representation<2>::SetBaseRank)
        .def("GetId", &Representation<2>::GetId)
        .def("GetLocalGrid", &Representation<2>::GetLocalGrid)
        .def("GetLocalWeights", &Representation<2>::GetLocalWeights)
        .def("MultiplyIntegrationWeights", Representation_2___MultiplyIntegrationWeightsWavefunction_2___Wavefunction_2___int)
        .def("MultiplyIntegrationWeights", Representation_2___MultiplyIntegrationWeightsWavefunction_2___Wavefunction_2__)
        .def("MultiplyIntegrationWeights", Representation_2___MultiplyIntegrationWeightsWavefunction_2___int)
        .def("MultiplyIntegrationWeights", Representation_2___MultiplyIntegrationWeightsWavefunction_2__)
        .def("GetGlobalOverlapMatrix", &Representation<2>::GetGlobalOverlapMatrix)
        .def("MultiplyOverlap", Representation_2___MultiplyOverlapWavefunction_2___Wavefunction_2___int)
        .def("MultiplyOverlap", Representation_2___MultiplyOverlapstd__complex_double__Wavefunction_2___std__complex_double__Wavefunction_2___int)
        .def("MultiplyOverlap", Representation_2___MultiplyOverlapWavefunction_2__)
        .def("SolveOverlap", &Representation<2>::SolveOverlap)
        .def("MultiplySqrtOverlap", &Representation<2>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &Representation<2>::SolveSqrtOverlap)
        .def("IsOrthogonalBasis", &Representation<2>::IsOrthogonalBasis)
        .def("GetFullShape", &Representation<2>::GetFullShape)
        .def("InnerProduct", &Representation<2>::InnerProduct)
        .def("GetGlobalWeights", &Representation<2>::GetGlobalWeights)
        .def("GetGlobalGrid", &Representation<2>::GetGlobalGrid)
        .def("ApplyConfigSection", &Representation<2>::ApplyConfigSection)
        .def("Copy", &Representation<2>::Copy)
    );
    register_ptr_to_python< boost::shared_ptr< Representation<2> > >();
    delete Representation_2_scope;

    scope* Representation_3_scope = new scope(
    class_< Representation<3>, boost::noncopyable >("Representation_3", no_init)
        .def("SetDistributedModel", &Representation<3>::SetDistributedModel)
        .def("GetDistributedModel", &Representation<3>::GetDistributedModel)
        .def("GetInitialShape", &Representation<3>::GetInitialShape)
        .def("GetBaseRank", &Representation<3>::GetBaseRank)
        .def("SetBaseRank", &Representation<3>::SetBaseRank)
        .def("GetId", &Representation<3>::GetId)
        .def("GetLocalGrid", &Representation<3>::GetLocalGrid)
        .def("GetLocalWeights", &Representation<3>::GetLocalWeights)
        .def("MultiplyIntegrationWeights", Representation_3___MultiplyIntegrationWeightsWavefunction_3___Wavefunction_3___int)
        .def("MultiplyIntegrationWeights", Representation_3___MultiplyIntegrationWeightsWavefunction_3___Wavefunction_3__)
        .def("MultiplyIntegrationWeights", Representation_3___MultiplyIntegrationWeightsWavefunction_3___int)
        .def("MultiplyIntegrationWeights", Representation_3___MultiplyIntegrationWeightsWavefunction_3__)
        .def("GetGlobalOverlapMatrix", &Representation<3>::GetGlobalOverlapMatrix)
        .def("MultiplyOverlap", Representation_3___MultiplyOverlapWavefunction_3___Wavefunction_3___int)
        .def("MultiplyOverlap", Representation_3___MultiplyOverlapstd__complex_double__Wavefunction_3___std__complex_double__Wavefunction_3___int)
        .def("MultiplyOverlap", Representation_3___MultiplyOverlapWavefunction_3__)
        .def("SolveOverlap", &Representation<3>::SolveOverlap)
        .def("MultiplySqrtOverlap", &Representation<3>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &Representation<3>::SolveSqrtOverlap)
        .def("IsOrthogonalBasis", &Representation<3>::IsOrthogonalBasis)
        .def("GetFullShape", &Representation<3>::GetFullShape)
        .def("InnerProduct", &Representation<3>::InnerProduct)
        .def("GetGlobalWeights", &Representation<3>::GetGlobalWeights)
        .def("GetGlobalGrid", &Representation<3>::GetGlobalGrid)
        .def("ApplyConfigSection", &Representation<3>::ApplyConfigSection)
        .def("Copy", &Representation<3>::Copy)
    );
    register_ptr_to_python< boost::shared_ptr< Representation<3> > >();
    delete Representation_3_scope;

    scope* Representation_4_scope = new scope(
    class_< Representation<4>, boost::noncopyable >("Representation_4", no_init)
        .def("SetDistributedModel", &Representation<4>::SetDistributedModel)
        .def("GetDistributedModel", &Representation<4>::GetDistributedModel)
        .def("GetInitialShape", &Representation<4>::GetInitialShape)
        .def("GetBaseRank", &Representation<4>::GetBaseRank)
        .def("SetBaseRank", &Representation<4>::SetBaseRank)
        .def("GetId", &Representation<4>::GetId)
        .def("GetLocalGrid", &Representation<4>::GetLocalGrid)
        .def("GetLocalWeights", &Representation<4>::GetLocalWeights)
        .def("MultiplyIntegrationWeights", Representation_4___MultiplyIntegrationWeightsWavefunction_4___Wavefunction_4___int)
        .def("MultiplyIntegrationWeights", Representation_4___MultiplyIntegrationWeightsWavefunction_4___Wavefunction_4__)
        .def("MultiplyIntegrationWeights", Representation_4___MultiplyIntegrationWeightsWavefunction_4___int)
        .def("MultiplyIntegrationWeights", Representation_4___MultiplyIntegrationWeightsWavefunction_4__)
        .def("GetGlobalOverlapMatrix", &Representation<4>::GetGlobalOverlapMatrix)
        .def("MultiplyOverlap", Representation_4___MultiplyOverlapWavefunction_4___Wavefunction_4___int)
        .def("MultiplyOverlap", Representation_4___MultiplyOverlapstd__complex_double__Wavefunction_4___std__complex_double__Wavefunction_4___int)
        .def("MultiplyOverlap", Representation_4___MultiplyOverlapWavefunction_4__)
        .def("SolveOverlap", &Representation<4>::SolveOverlap)
        .def("MultiplySqrtOverlap", &Representation<4>::MultiplySqrtOverlap)
        .def("SolveSqrtOverlap", &Representation<4>::SolveSqrtOverlap)
        .def("IsOrthogonalBasis", &Representation<4>::IsOrthogonalBasis)
        .def("GetFullShape", &Representation<4>::GetFullShape)
        .def("InnerProduct", &Representation<4>::InnerProduct)
        .def("GetGlobalWeights", &Representation<4>::GetGlobalWeights)
        .def("GetGlobalGrid", &Representation<4>::GetGlobalGrid)
        .def("ApplyConfigSection", &Representation<4>::ApplyConfigSection)
        .def("Copy", &Representation<4>::Copy)
    );
    register_ptr_to_python< boost::shared_ptr< Representation<4> > >();
    delete Representation_4_scope;

}


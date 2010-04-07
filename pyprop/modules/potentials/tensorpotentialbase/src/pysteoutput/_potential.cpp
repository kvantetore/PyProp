
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <src/potential.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_pyprop_modules_potentials_tensorpotentialbase_src_potential()
{
    class_< DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>, boost::noncopyable >("DynamicPotentialEvaluator_KineticEnergyPotential_1_1", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>, boost::noncopyable >("DynamicPotentialEvaluator_KineticEnergyPotential_2_2", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>, boost::noncopyable >("DynamicPotentialEvaluator_KineticEnergyPotential_3_3", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<KineticEnergyPotential<3>,3>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<KineticEnergyPotential<4>,4>, boost::noncopyable >("DynamicPotentialEvaluator_KineticEnergyPotential_4_4", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<KineticEnergyPotential<4>,4>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<4>,4>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<4>,4>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<4>,4>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<4>,4>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<KineticEnergyPotential<4>,4>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<KineticEnergyPotential<4>,4>::CalculateExpectationValue)
    ;

}


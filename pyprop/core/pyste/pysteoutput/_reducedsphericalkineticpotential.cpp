
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential/reducedsphericalkineticpotential.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_reducedsphericalkineticpotential()
{
    class_< DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<2>,2>, boost::noncopyable >("ReducedAngularKineticEnergyPotential_2", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<3>,3>, boost::noncopyable >("ReducedAngularKineticEnergyPotential_3", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<ReducedAngularKineticEnergyPotential<3>,3>::CalculateExpectationValue)
    ;

}


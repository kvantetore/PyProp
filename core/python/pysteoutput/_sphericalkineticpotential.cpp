
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential/sphericalkineticenergypotential.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_sphericalkineticpotential()
{
    class_< SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<3>,2> >("AngularKineticEnergyPotential_2", init<  >())
        .def(init< const SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<3>,2>& >())
        .def("ApplyConfigSection", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<3>,2>::ApplyConfigSection)
        .def("ApplyPotential", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<3>,2>::ApplyPotential)
        .def("MultiplyPotential", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<3>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<3>,2>::UpdateStaticPotential)
        .def("GetPotential", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<3>,2>::GetPotential)
        .def("CalculateExpectationValue", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<3>,2>::CalculateExpectationValue)
    ;

    class_< SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<4>,3> >("AngularKineticEnergyPotential_3", init<  >())
        .def(init< const SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<4>,3>& >())
        .def("ApplyConfigSection", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<4>,3>::ApplyConfigSection)
        .def("ApplyPotential", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<4>,3>::ApplyPotential)
        .def("MultiplyPotential", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<4>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<4>,3>::UpdateStaticPotential)
        .def("GetPotential", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<4>,3>::GetPotential)
        .def("CalculateExpectationValue", &SphericalDynamicPotentialEvaluator<AngularKineticEnergyPotential<4>,3>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<RadialKineticEnergyPotential<1>,1> >("RadialKineticEnergyPotential_1", init<  >())
        .def(init< const DynamicPotentialEvaluator<RadialKineticEnergyPotential<1>,1>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<1>,1>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<1>,1>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<1>,1>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<1>,1>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<1>,1>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<RadialKineticEnergyPotential<2>,2> >("RadialKineticEnergyPotential_2", init<  >())
        .def(init< const DynamicPotentialEvaluator<RadialKineticEnergyPotential<2>,2>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<RadialKineticEnergyPotential<3>,3> >("RadialKineticEnergyPotential_3", init<  >())
        .def(init< const DynamicPotentialEvaluator<RadialKineticEnergyPotential<3>,3>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<RadialKineticEnergyPotential<3>,3>::CalculateExpectationValue)
    ;

}


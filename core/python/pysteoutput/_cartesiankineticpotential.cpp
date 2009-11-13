
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential/cartesiankineticenergypotential.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_cartesiankineticpotential()
{
    class_< DynamicPotentialEvaluator<CartesianKineticEnergyPotential<1>,1>, boost::noncopyable >("DynamicPotentialEvaluator_CartesianKineticEnergyPotential_1_1", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<1>,1>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<1>,1>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<1>,1>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<1>,1>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<1>,1>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<CartesianKineticEnergyPotential<2>,2>, boost::noncopyable >("DynamicPotentialEvaluator_CartesianKineticEnergyPotential_2_2", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<CartesianKineticEnergyPotential<3>,3>, boost::noncopyable >("DynamicPotentialEvaluator_CartesianKineticEnergyPotential_3_3", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<3>,3>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<CartesianKineticEnergyPotential<4>,4>, boost::noncopyable >("DynamicPotentialEvaluator_CartesianKineticEnergyPotential_4_4", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<4>,4>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<4>,4>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<4>,4>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<4>,4>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<4>,4>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<4>,4>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<CartesianKineticEnergyPotential<4>,4>::CalculateExpectationValue)
    ;

}


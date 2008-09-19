
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential/polarpotential.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_polarpotential()
{
    class_< DynamicPotentialEvaluator<PolarKineticPotential<2>,2> >("PolarKineticPotential_2", init<  >())
        .def(init< const DynamicPotentialEvaluator<PolarKineticPotential<2>,2>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<PolarKineticPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<PolarKineticPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<PolarKineticPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<PolarKineticPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<PolarKineticPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<PolarKineticPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<PolarKineticPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<PolarKineticPotential<3>,3> >("PolarKineticPotential_3", init<  >())
        .def(init< const DynamicPotentialEvaluator<PolarKineticPotential<3>,3>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<PolarKineticPotential<3>,3>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<PolarKineticPotential<3>,3>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<PolarKineticPotential<3>,3>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<PolarKineticPotential<3>,3>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<PolarKineticPotential<3>,3>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<PolarKineticPotential<3>,3>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<PolarKineticPotential<3>,3>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<PolarKineticPotential<4>,4> >("PolarKineticPotential_4", init<  >())
        .def(init< const DynamicPotentialEvaluator<PolarKineticPotential<4>,4>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<PolarKineticPotential<4>,4>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<PolarKineticPotential<4>,4>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<PolarKineticPotential<4>,4>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<PolarKineticPotential<4>,4>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<PolarKineticPotential<4>,4>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<PolarKineticPotential<4>,4>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<PolarKineticPotential<4>,4>::CalculateExpectationValue)
    ;

}


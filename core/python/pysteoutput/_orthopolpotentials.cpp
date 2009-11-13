
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential/orthopolpotentials.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_orthopolpotentials()
{
    class_< RankOnePotentialEvaluator<AddedHermitePotential<1>,1>, boost::noncopyable >("RankOnePotentialEvaluator_AddedHermitePotential_1_1", no_init)
        .def("ApplyConfigSection", &RankOnePotentialEvaluator<AddedHermitePotential<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,1>::MultiplyPotential)
        .def("GetPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,1>::GetPotential)
        .def("CalculateExpectationValue", &RankOnePotentialEvaluator<AddedHermitePotential<1>,1>::CalculateExpectationValue)
    ;

    class_< RankOnePotentialEvaluator<AddedHermitePotential<1>,2>, boost::noncopyable >("RankOnePotentialEvaluator_AddedHermitePotential_1_2", no_init)
        .def("ApplyConfigSection", &RankOnePotentialEvaluator<AddedHermitePotential<1>,2>::ApplyConfigSection)
        .def("ApplyPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,2>::ApplyPotential)
        .def("MultiplyPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,2>::MultiplyPotential)
        .def("GetPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,2>::GetPotential)
        .def("CalculateExpectationValue", &RankOnePotentialEvaluator<AddedHermitePotential<1>,2>::CalculateExpectationValue)
    ;

    class_< RankOnePotentialEvaluator<AddedHermitePotential<1>,3>, boost::noncopyable >("RankOnePotentialEvaluator_AddedHermitePotential_1_3", no_init)
        .def("ApplyConfigSection", &RankOnePotentialEvaluator<AddedHermitePotential<1>,3>::ApplyConfigSection)
        .def("ApplyPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,3>::ApplyPotential)
        .def("MultiplyPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,3>::MultiplyPotential)
        .def("GetPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,3>::GetPotential)
        .def("CalculateExpectationValue", &RankOnePotentialEvaluator<AddedHermitePotential<1>,3>::CalculateExpectationValue)
    ;

    class_< RankOnePotentialEvaluator<AddedHermitePotential<1>,4>, boost::noncopyable >("RankOnePotentialEvaluator_AddedHermitePotential_1_4", no_init)
        .def("ApplyConfigSection", &RankOnePotentialEvaluator<AddedHermitePotential<1>,4>::ApplyConfigSection)
        .def("ApplyPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,4>::ApplyPotential)
        .def("MultiplyPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,4>::MultiplyPotential)
        .def("GetPotential", &RankOnePotentialEvaluator<AddedHermitePotential<1>,4>::GetPotential)
        .def("CalculateExpectationValue", &RankOnePotentialEvaluator<AddedHermitePotential<1>,4>::CalculateExpectationValue)
    ;

    class_< RankOnePotentialEvaluator<AddedLaguerrePotential<1>,1>, boost::noncopyable >("RankOnePotentialEvaluator_AddedLaguerrePotential_1_1", no_init)
        .def("ApplyConfigSection", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,1>::MultiplyPotential)
        .def("GetPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,1>::GetPotential)
        .def("CalculateExpectationValue", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,1>::CalculateExpectationValue)
    ;

    class_< RankOnePotentialEvaluator<AddedLaguerrePotential<1>,2>, boost::noncopyable >("RankOnePotentialEvaluator_AddedLaguerrePotential_1_2", no_init)
        .def("ApplyConfigSection", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,2>::ApplyConfigSection)
        .def("ApplyPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,2>::ApplyPotential)
        .def("MultiplyPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,2>::MultiplyPotential)
        .def("GetPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,2>::GetPotential)
        .def("CalculateExpectationValue", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,2>::CalculateExpectationValue)
    ;

    class_< RankOnePotentialEvaluator<AddedLaguerrePotential<1>,3>, boost::noncopyable >("RankOnePotentialEvaluator_AddedLaguerrePotential_1_3", no_init)
        .def("ApplyConfigSection", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,3>::ApplyConfigSection)
        .def("ApplyPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,3>::ApplyPotential)
        .def("MultiplyPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,3>::MultiplyPotential)
        .def("GetPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,3>::GetPotential)
        .def("CalculateExpectationValue", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,3>::CalculateExpectationValue)
    ;

    class_< RankOnePotentialEvaluator<AddedLaguerrePotential<1>,4>, boost::noncopyable >("RankOnePotentialEvaluator_AddedLaguerrePotential_1_4", no_init)
        .def("ApplyConfigSection", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,4>::ApplyConfigSection)
        .def("ApplyPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,4>::ApplyPotential)
        .def("MultiplyPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,4>::MultiplyPotential)
        .def("GetPotential", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,4>::GetPotential)
        .def("CalculateExpectationValue", &RankOnePotentialEvaluator<AddedLaguerrePotential<1>,4>::CalculateExpectationValue)
    ;

}


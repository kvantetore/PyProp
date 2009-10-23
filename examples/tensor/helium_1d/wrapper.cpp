
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpotential)
{
    class_< DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>, boost::noncopyable >("KineticEnergyPotential_1", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<KineticEnergyPotential<1>,1>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>, boost::noncopyable >("KineticEnergyPotential_2", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<KineticEnergyPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<DipoleLaserPotentialLength<2>,2>, boost::noncopyable >("DipoleLaserPotentialLength_2", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<DipoleLaserPotentialLength<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<DipoleLaserPotentialLength<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<DipoleLaserPotentialLength<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<DipoleLaserPotentialLength<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<DipoleLaserPotentialLength<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<DipoleLaserPotentialLength<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<DipoleLaserPotentialLength<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<SoftCoulombPotential<2>,2>, boost::noncopyable >("SoftCoulombPotential_2", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<SoftCoulombPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<SoftCoulombPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<SoftCoulombPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<SoftCoulombPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<SoftCoulombPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<SoftCoulombPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<SoftCoulombPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<SoftCoulombPotential1D<1>,1>, boost::noncopyable >("SoftCoulombPotential1D_1", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<SoftCoulombPotential1D<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<SoftCoulombPotential1D<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<SoftCoulombPotential1D<1>,1>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<SoftCoulombPotential1D<1>,1>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<SoftCoulombPotential1D<1>,1>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<SoftCoulombPotential1D<1>,1>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<SoftCoulombPotential1D<1>,1>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<TwoElectronCorrelation1D<2>,2>, boost::noncopyable >("TwoElectronCorrelation1D_2", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<TwoElectronCorrelation1D<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<TwoElectronCorrelation1D<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<TwoElectronCorrelation1D<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<TwoElectronCorrelation1D<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<TwoElectronCorrelation1D<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<TwoElectronCorrelation1D<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<TwoElectronCorrelation1D<2>,2>::CalculateExpectationValue)
    ;

}


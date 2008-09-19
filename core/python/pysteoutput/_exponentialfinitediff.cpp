
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <finitediff/exponentialfinitedifference.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_exponentialfinitediff()
{
    def("EXPFD_GetDiagonalValue", &EXPFD_GetDiagonalValue);
    def("EXPFD_UpdateEigenvalues", &EXPFD_UpdateEigenvalues);
    def("EXPFD_UpdateEigenvectors", &EXPFD_UpdateEigenvectors);
    def("EXPFD_UpdateBlock", &EXPFD_UpdateBlock);
    def("EXPFD_UpdateDiagonal", &EXPFD_UpdateDiagonal);
    class_< ExponentialFiniteDifferenceEvaluator<HarmonicOscillatorPotential<1>,1> >("ExponentialFiniteDifferenceEvaluator_HarmonicOscillatorPotential_1_1", init<  >())
        .def(init< const ExponentialFiniteDifferenceEvaluator<HarmonicOscillatorPotential<1>,1>& >())
        .def_readwrite("strength", &HarmonicOscillatorPotential<1>::strength)
        .def_readwrite("CurTime", &PotentialBase<1>::CurTime)
        .def_readwrite("TimeStep", &PotentialBase<1>::TimeStep)
        .def("UpdateWavefunction", &ExponentialFiniteDifferenceEvaluator<HarmonicOscillatorPotential<1>,1>::UpdateWavefunction)
        .def("MultiplyOperator", &ExponentialFiniteDifferenceEvaluator<HarmonicOscillatorPotential<1>,1>::MultiplyOperator)
        .def("ApplyConfigSection", &HarmonicOscillatorPotential<1>::ApplyConfigSection)
        .def("GetPotentialValue", &HarmonicOscillatorPotential<1>::GetPotentialValue)
        .def("CurTimeUpdated", &PotentialBase<1>::CurTimeUpdated)
        .def("IsTimeDependent", &PotentialBase<1>::IsTimeDependent)
    ;

}


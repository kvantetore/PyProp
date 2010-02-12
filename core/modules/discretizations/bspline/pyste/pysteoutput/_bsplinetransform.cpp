
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <combinedpropagator/bsplinepropagator.h>
#include <combinedpropagator/bsplinetransform.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

void (BSpline::Propagator<1>::*BSpline__Propagator_1___Setupconststd__complex_double___constWavefunction_1___boost__shared_ptr_BSpline__BSpline__int)(const std::complex<double>&, const Wavefunction<1>&, boost::shared_ptr<BSpline::BSpline>, int)  = &BSpline::Propagator<1>::Setup;

void (BSpline::Propagator<1>::*BSpline__Propagator_1___Setupconststd__complex_double___constWavefunction_1___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)(const std::complex<double>&, const Wavefunction<1>&, boost::shared_ptr<BSpline::BSpline>, blitz::Array<double,1>, int)  = &BSpline::Propagator<1>::Setup;

void (BSpline::Propagator<2>::*BSpline__Propagator_2___Setupconststd__complex_double___constWavefunction_2___boost__shared_ptr_BSpline__BSpline__int)(const std::complex<double>&, const Wavefunction<2>&, boost::shared_ptr<BSpline::BSpline>, int)  = &BSpline::Propagator<2>::Setup;

void (BSpline::Propagator<2>::*BSpline__Propagator_2___Setupconststd__complex_double___constWavefunction_2___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)(const std::complex<double>&, const Wavefunction<2>&, boost::shared_ptr<BSpline::BSpline>, blitz::Array<double,1>, int)  = &BSpline::Propagator<2>::Setup;

void (BSpline::Propagator<3>::*BSpline__Propagator_3___Setupconststd__complex_double___constWavefunction_3___boost__shared_ptr_BSpline__BSpline__int)(const std::complex<double>&, const Wavefunction<3>&, boost::shared_ptr<BSpline::BSpline>, int)  = &BSpline::Propagator<3>::Setup;

void (BSpline::Propagator<3>::*BSpline__Propagator_3___Setupconststd__complex_double___constWavefunction_3___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)(const std::complex<double>&, const Wavefunction<3>&, boost::shared_ptr<BSpline::BSpline>, blitz::Array<double,1>, int)  = &BSpline::Propagator<3>::Setup;

void (BSpline::Propagator<4>::*BSpline__Propagator_4___Setupconststd__complex_double___constWavefunction_4___boost__shared_ptr_BSpline__BSpline__int)(const std::complex<double>&, const Wavefunction<4>&, boost::shared_ptr<BSpline::BSpline>, int)  = &BSpline::Propagator<4>::Setup;

void (BSpline::Propagator<4>::*BSpline__Propagator_4___Setupconststd__complex_double___constWavefunction_4___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)(const std::complex<double>&, const Wavefunction<4>&, boost::shared_ptr<BSpline::BSpline>, blitz::Array<double,1>, int)  = &BSpline::Propagator<4>::Setup;


}// namespace 


// Module ======================================================================
void Export_core_modules_discretizations_bspline_pyste_bsplinetransform()
{
    class_< BSpline::Propagator<1>, boost::noncopyable >("BSplinePropagator_1", init<  >())
        .def("ApplyConfigSection", &BSpline::Propagator<1>::ApplyConfigSection)
        .def("Setup", BSpline__Propagator_1___Setupconststd__complex_double___constWavefunction_1___boost__shared_ptr_BSpline__BSpline__int)
        .def("Setup", BSpline__Propagator_1___Setupconststd__complex_double___constWavefunction_1___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)
        .def("SetupCentrifugalPotential", &BSpline::Propagator<1>::SetupCentrifugalPotential)
        .def("AdvanceStep", &BSpline::Propagator<1>::AdvanceStep)
        .def("MultiplyHamiltonian", &BSpline::Propagator<1>::MultiplyHamiltonian)
        .def("ApplyCrankNicolson", &BSpline::Propagator<1>::ApplyCrankNicolson)
        .def("GetPotentialSlice", &BSpline::Propagator<1>::GetPotentialSlice)
        .def("GetOverlapMatrix", &BSpline::Propagator<1>::GetOverlapMatrix)
        .def("GetHamiltonianMatrix", &BSpline::Propagator<1>::GetHamiltonianMatrix)
        .def("GetOverlapMatrixBlas", &BSpline::Propagator<1>::GetOverlapMatrixBlas)
        .def("GetPropagationMatrix", &BSpline::Propagator<1>::GetPropagationMatrix)
        .def("GetCentrifugalMatrixBlas", &BSpline::Propagator<1>::GetCentrifugalMatrixBlas)
        .def("GetCentrifugalMatrix", &BSpline::Propagator<1>::GetCentrifugalMatrix)
        .def("GetBigPropagationMatrix", &BSpline::Propagator<1>::GetBigPropagationMatrix)
        .def("GetGlobalLmax", &BSpline::Propagator<1>::GetGlobalLmax)
        .def("SetPropagationAlgorithm", &BSpline::Propagator<1>::SetPropagationAlgorithm)
        .def("GetPropagationAlgorithm", &BSpline::Propagator<1>::GetPropagationAlgorithm)
    ;

    class_< BSpline::Propagator<2>, boost::noncopyable >("BSplinePropagator_2", init<  >())
        .def("ApplyConfigSection", &BSpline::Propagator<2>::ApplyConfigSection)
        .def("Setup", BSpline__Propagator_2___Setupconststd__complex_double___constWavefunction_2___boost__shared_ptr_BSpline__BSpline__int)
        .def("Setup", BSpline__Propagator_2___Setupconststd__complex_double___constWavefunction_2___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)
        .def("SetupCentrifugalPotential", &BSpline::Propagator<2>::SetupCentrifugalPotential)
        .def("AdvanceStep", &BSpline::Propagator<2>::AdvanceStep)
        .def("MultiplyHamiltonian", &BSpline::Propagator<2>::MultiplyHamiltonian)
        .def("ApplyCrankNicolson", &BSpline::Propagator<2>::ApplyCrankNicolson)
        .def("GetPotentialSlice", &BSpline::Propagator<2>::GetPotentialSlice)
        .def("GetOverlapMatrix", &BSpline::Propagator<2>::GetOverlapMatrix)
        .def("GetHamiltonianMatrix", &BSpline::Propagator<2>::GetHamiltonianMatrix)
        .def("GetOverlapMatrixBlas", &BSpline::Propagator<2>::GetOverlapMatrixBlas)
        .def("GetPropagationMatrix", &BSpline::Propagator<2>::GetPropagationMatrix)
        .def("GetCentrifugalMatrixBlas", &BSpline::Propagator<2>::GetCentrifugalMatrixBlas)
        .def("GetCentrifugalMatrix", &BSpline::Propagator<2>::GetCentrifugalMatrix)
        .def("GetBigPropagationMatrix", &BSpline::Propagator<2>::GetBigPropagationMatrix)
        .def("GetGlobalLmax", &BSpline::Propagator<2>::GetGlobalLmax)
        .def("SetPropagationAlgorithm", &BSpline::Propagator<2>::SetPropagationAlgorithm)
        .def("GetPropagationAlgorithm", &BSpline::Propagator<2>::GetPropagationAlgorithm)
    ;

    class_< BSpline::Propagator<3>, boost::noncopyable >("BSplinePropagator_3", init<  >())
        .def("ApplyConfigSection", &BSpline::Propagator<3>::ApplyConfigSection)
        .def("Setup", BSpline__Propagator_3___Setupconststd__complex_double___constWavefunction_3___boost__shared_ptr_BSpline__BSpline__int)
        .def("Setup", BSpline__Propagator_3___Setupconststd__complex_double___constWavefunction_3___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)
        .def("SetupCentrifugalPotential", &BSpline::Propagator<3>::SetupCentrifugalPotential)
        .def("AdvanceStep", &BSpline::Propagator<3>::AdvanceStep)
        .def("MultiplyHamiltonian", &BSpline::Propagator<3>::MultiplyHamiltonian)
        .def("ApplyCrankNicolson", &BSpline::Propagator<3>::ApplyCrankNicolson)
        .def("GetPotentialSlice", &BSpline::Propagator<3>::GetPotentialSlice)
        .def("GetOverlapMatrix", &BSpline::Propagator<3>::GetOverlapMatrix)
        .def("GetHamiltonianMatrix", &BSpline::Propagator<3>::GetHamiltonianMatrix)
        .def("GetOverlapMatrixBlas", &BSpline::Propagator<3>::GetOverlapMatrixBlas)
        .def("GetPropagationMatrix", &BSpline::Propagator<3>::GetPropagationMatrix)
        .def("GetCentrifugalMatrixBlas", &BSpline::Propagator<3>::GetCentrifugalMatrixBlas)
        .def("GetCentrifugalMatrix", &BSpline::Propagator<3>::GetCentrifugalMatrix)
        .def("GetBigPropagationMatrix", &BSpline::Propagator<3>::GetBigPropagationMatrix)
        .def("GetGlobalLmax", &BSpline::Propagator<3>::GetGlobalLmax)
        .def("SetPropagationAlgorithm", &BSpline::Propagator<3>::SetPropagationAlgorithm)
        .def("GetPropagationAlgorithm", &BSpline::Propagator<3>::GetPropagationAlgorithm)
    ;

    class_< BSpline::Propagator<4>, boost::noncopyable >("BSplinePropagator_4", init<  >())
        .def("ApplyConfigSection", &BSpline::Propagator<4>::ApplyConfigSection)
        .def("Setup", BSpline__Propagator_4___Setupconststd__complex_double___constWavefunction_4___boost__shared_ptr_BSpline__BSpline__int)
        .def("Setup", BSpline__Propagator_4___Setupconststd__complex_double___constWavefunction_4___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)
        .def("SetupCentrifugalPotential", &BSpline::Propagator<4>::SetupCentrifugalPotential)
        .def("AdvanceStep", &BSpline::Propagator<4>::AdvanceStep)
        .def("MultiplyHamiltonian", &BSpline::Propagator<4>::MultiplyHamiltonian)
        .def("ApplyCrankNicolson", &BSpline::Propagator<4>::ApplyCrankNicolson)
        .def("GetPotentialSlice", &BSpline::Propagator<4>::GetPotentialSlice)
        .def("GetOverlapMatrix", &BSpline::Propagator<4>::GetOverlapMatrix)
        .def("GetHamiltonianMatrix", &BSpline::Propagator<4>::GetHamiltonianMatrix)
        .def("GetOverlapMatrixBlas", &BSpline::Propagator<4>::GetOverlapMatrixBlas)
        .def("GetPropagationMatrix", &BSpline::Propagator<4>::GetPropagationMatrix)
        .def("GetCentrifugalMatrixBlas", &BSpline::Propagator<4>::GetCentrifugalMatrixBlas)
        .def("GetCentrifugalMatrix", &BSpline::Propagator<4>::GetCentrifugalMatrix)
        .def("GetBigPropagationMatrix", &BSpline::Propagator<4>::GetBigPropagationMatrix)
        .def("GetGlobalLmax", &BSpline::Propagator<4>::GetGlobalLmax)
        .def("SetPropagationAlgorithm", &BSpline::Propagator<4>::SetPropagationAlgorithm)
        .def("GetPropagationAlgorithm", &BSpline::Propagator<4>::GetPropagationAlgorithm)
    ;

    class_< BSpline::BSplineTransform<1>, boost::noncopyable >("BSplineTransform_1", init<  >())
        .def("SetupStep", &BSpline::BSplineTransform<1>::SetupStep)
        .def("ForwardTransform", &BSpline::BSplineTransform<1>::ForwardTransform)
        .def("InverseTransform", &BSpline::BSplineTransform<1>::InverseTransform)
        .def("GetBaseRank", &BSpline::BSplineTransform<1>::GetBaseRank)
        .def("SetBaseRank", &BSpline::BSplineTransform<1>::SetBaseRank)
    ;

    class_< BSpline::BSplineTransform<2>, boost::noncopyable >("BSplineTransform_2", init<  >())
        .def("SetupStep", &BSpline::BSplineTransform<2>::SetupStep)
        .def("ForwardTransform", &BSpline::BSplineTransform<2>::ForwardTransform)
        .def("InverseTransform", &BSpline::BSplineTransform<2>::InverseTransform)
        .def("GetBaseRank", &BSpline::BSplineTransform<2>::GetBaseRank)
        .def("SetBaseRank", &BSpline::BSplineTransform<2>::SetBaseRank)
    ;

    class_< BSpline::BSplineTransform<3>, boost::noncopyable >("BSplineTransform_3", init<  >())
        .def("SetupStep", &BSpline::BSplineTransform<3>::SetupStep)
        .def("ForwardTransform", &BSpline::BSplineTransform<3>::ForwardTransform)
        .def("InverseTransform", &BSpline::BSplineTransform<3>::InverseTransform)
        .def("GetBaseRank", &BSpline::BSplineTransform<3>::GetBaseRank)
        .def("SetBaseRank", &BSpline::BSplineTransform<3>::SetBaseRank)
    ;

    class_< BSpline::BSplineTransform<4>, boost::noncopyable >("BSplineTransform_4", init<  >())
        .def("SetupStep", &BSpline::BSplineTransform<4>::SetupStep)
        .def("ForwardTransform", &BSpline::BSplineTransform<4>::ForwardTransform)
        .def("InverseTransform", &BSpline::BSplineTransform<4>::InverseTransform)
        .def("GetBaseRank", &BSpline::BSplineTransform<4>::GetBaseRank)
        .def("SetBaseRank", &BSpline::BSplineTransform<4>::SetBaseRank)
    ;

}


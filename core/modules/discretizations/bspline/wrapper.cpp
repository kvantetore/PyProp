
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <../representation/bspline/bsplinegridrepresentation.h>
#include <../representation/bspline/bsplinerepresentation.h>
#include <../transform/bspline/bspline.h>
#include <../transform/bspline/bsplinepropagator.h>
#include <../transform/bspline/bsplinesolver.h>
#include <../transform/bspline/bsplinetransform.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
#include "../tensorpotential/basis_bspline.h"


namespace  {

void (BSpline::Propagator<1>::*BSpline__Propagator_1___Setupconststd__complex_double___constWavefunction_1___boost__shared_ptr_BSpline__BSpline__int)(const std::complex<double>&, const Wavefunction<1>&, boost::shared_ptr<BSpline::BSpline>, int)  = &BSpline::Propagator<1>::Setup;

void (BSpline::Propagator<1>::*BSpline__Propagator_1___Setupconststd__complex_double___constWavefunction_1___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)(const std::complex<double>&, const Wavefunction<1>&, boost::shared_ptr<BSpline::BSpline>, blitz::Array<double,1>, int)  = &BSpline::Propagator<1>::Setup;

void (BSpline::Propagator<2>::*BSpline__Propagator_2___Setupconststd__complex_double___constWavefunction_2___boost__shared_ptr_BSpline__BSpline__int)(const std::complex<double>&, const Wavefunction<2>&, boost::shared_ptr<BSpline::BSpline>, int)  = &BSpline::Propagator<2>::Setup;

void (BSpline::Propagator<2>::*BSpline__Propagator_2___Setupconststd__complex_double___constWavefunction_2___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)(const std::complex<double>&, const Wavefunction<2>&, boost::shared_ptr<BSpline::BSpline>, blitz::Array<double,1>, int)  = &BSpline::Propagator<2>::Setup;

void (BSpline::Propagator<3>::*BSpline__Propagator_3___Setupconststd__complex_double___constWavefunction_3___boost__shared_ptr_BSpline__BSpline__int)(const std::complex<double>&, const Wavefunction<3>&, boost::shared_ptr<BSpline::BSpline>, int)  = &BSpline::Propagator<3>::Setup;

void (BSpline::Propagator<3>::*BSpline__Propagator_3___Setupconststd__complex_double___constWavefunction_3___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)(const std::complex<double>&, const Wavefunction<3>&, boost::shared_ptr<BSpline::BSpline>, blitz::Array<double,1>, int)  = &BSpline::Propagator<3>::Setup;

void (BSpline::Propagator<4>::*BSpline__Propagator_4___Setupconststd__complex_double___constWavefunction_4___boost__shared_ptr_BSpline__BSpline__int)(const std::complex<double>&, const Wavefunction<4>&, boost::shared_ptr<BSpline::BSpline>, int)  = &BSpline::Propagator<4>::Setup;

void (BSpline::Propagator<4>::*BSpline__Propagator_4___Setupconststd__complex_double___constWavefunction_4___boost__shared_ptr_BSpline__BSpline__blitz__Array_double_1__int)(const std::complex<double>&, const Wavefunction<4>&, boost::shared_ptr<BSpline::BSpline>, blitz::Array<double,1>, int)  = &BSpline::Propagator<4>::Setup;

double (BSpline::BSpline::*BSpline__BSpline__BSplineOverlapIntegralint_int)(int, int)  = &BSpline::BSpline::BSplineOverlapIntegral;

double (BSpline::BSpline::*BSpline__BSpline__BSplineOverlapIntegralblitz__Array_double_1__int_int)(blitz::Array<double,1>, int, int)  = &BSpline::BSpline::BSplineOverlapIntegral;

blitz::Array<std::complex<double>,1> (BSpline::BSpline::*BSpline__BSpline__ExpandFunctionInBSplinesboost__python__object)(boost::python::object)  = &BSpline::BSpline::ExpandFunctionInBSplines;

void (BSpline::BSpline::*BSpline__BSpline__ExpandFunctionInBSplinesblitz__Array_std__complex_double__1__blitz__Array_std__complex_double__1_)(blitz::Array<std::complex<double>,1>, blitz::Array<std::complex<double>,1>)  = &BSpline::BSpline::ExpandFunctionInBSplines;

blitz::Array<std::complex<double>,1> (BSpline::BSpline::*BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1_)(blitz::Array<std::complex<double>,1>)  = &BSpline::BSpline::ConstructFunctionFromBSplineExpansion;

blitz::Array<std::complex<double>,1> (BSpline::BSpline::*BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1__blitz__Array_double_1_)(blitz::Array<std::complex<double>,1>, blitz::Array<double,1>)  = &BSpline::BSpline::ConstructFunctionFromBSplineExpansion;

void (BSpline::BSpline::*BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1__blitz__Array_double_1__blitz__Array_std__complex_double__1_)(blitz::Array<std::complex<double>,1>, blitz::Array<double,1>, blitz::Array<std::complex<double>,1>)  = &BSpline::BSpline::ConstructFunctionFromBSplineExpansion;


}// namespace 


// Module ======================================================================
BOOST_PYTHON_MODULE(libbspline)
{
    class_< BSpline::BSplineSolver<1>, boost::noncopyable >("BSpline_BSplineSolver_1", init<  >())
        .def_readwrite("MatrixData", &BSpline::BSplineSolver<1>::MatrixData)
        .def_readwrite("PotentialData", &BSpline::BSplineSolver<1>::PotentialData)
        .def_readwrite("PivotData", &BSpline::BSplineSolver<1>::PivotData)
        .def("AddTensorPotential", &BSpline::BSplineSolver<1>::AddTensorPotential)
        .def("Setup", &BSpline::BSplineSolver<1>::Setup)
        .def("Solve", &BSpline::BSplineSolver<1>::Solve)
    ;

    class_< BSpline::BSplineSolver<2>, boost::noncopyable >("BSpline_BSplineSolver_2", init<  >())
        .def_readwrite("MatrixData", &BSpline::BSplineSolver<2>::MatrixData)
        .def_readwrite("PotentialData", &BSpline::BSplineSolver<2>::PotentialData)
        .def_readwrite("PivotData", &BSpline::BSplineSolver<2>::PivotData)
        .def("AddTensorPotential", &BSpline::BSplineSolver<2>::AddTensorPotential)
        .def("Setup", &BSpline::BSplineSolver<2>::Setup)
        .def("Solve", &BSpline::BSplineSolver<2>::Solve)
    ;

    class_< BSpline::BSplineSolver<3>, boost::noncopyable >("BSpline_BSplineSolver_3", init<  >())
        .def_readwrite("MatrixData", &BSpline::BSplineSolver<3>::MatrixData)
        .def_readwrite("PotentialData", &BSpline::BSplineSolver<3>::PotentialData)
        .def_readwrite("PivotData", &BSpline::BSplineSolver<3>::PivotData)
        .def("AddTensorPotential", &BSpline::BSplineSolver<3>::AddTensorPotential)
        .def("Setup", &BSpline::BSplineSolver<3>::Setup)
        .def("Solve", &BSpline::BSplineSolver<3>::Solve)
    ;

    class_< BSpline::BSplineSolver<4>, boost::noncopyable >("BSpline_BSplineSolver_4", init<  >())
        .def_readwrite("MatrixData", &BSpline::BSplineSolver<4>::MatrixData)
        .def_readwrite("PotentialData", &BSpline::BSplineSolver<4>::PotentialData)
        .def_readwrite("PivotData", &BSpline::BSplineSolver<4>::PivotData)
        .def("AddTensorPotential", &BSpline::BSplineSolver<4>::AddTensorPotential)
        .def("Setup", &BSpline::BSplineSolver<4>::Setup)
        .def("Solve", &BSpline::BSplineSolver<4>::Solve)
    ;

    scope* BSpline_BSplineRepresentation_scope = new scope(
    class_< BSpline::BSplineRepresentation, bases< Representation<1> >  >("BSplineRepresentation", init<  >())
        .def(init< const BSpline::BSplineRepresentation& >())
        .def("Copy", &BSpline::BSplineRepresentation::Copy)
        .def("GetFullShape", &BSpline::BSplineRepresentation::GetFullShape)
        .def("InnerProduct", &BSpline::BSplineRepresentation::InnerProduct)
        .def("GetGlobalWeights", &BSpline::BSplineRepresentation::GetGlobalWeights)
        .def("GetGlobalGrid", &BSpline::BSplineRepresentation::GetGlobalGrid)
        .def("GetGlobalOverlapMatrix", &BSpline::BSplineRepresentation::GetGlobalOverlapMatrix)
        .def("ApplyConfigSection", &BSpline::BSplineRepresentation::ApplyConfigSection)
        .def("SetupRepresentation", &BSpline::BSplineRepresentation::SetupRepresentation)
        .def("GetBSplineObject", &BSpline::BSplineRepresentation::GetBSplineObject)
    );
    register_ptr_to_python< boost::shared_ptr< BSpline::BSplineRepresentation > >();
    delete BSpline_BSplineRepresentation_scope;

    scope* BSpline_BSplineGridRepresentation_scope = new scope(
    class_< BSpline::BSplineGridRepresentation, bases< Representation<1> >  >("BSplineGridRepresentation", init<  >())
        .def(init< const BSpline::BSplineGridRepresentation& >())
        .def("Copy", &BSpline::BSplineGridRepresentation::Copy)
        .def("GetFullShape", &BSpline::BSplineGridRepresentation::GetFullShape)
        .def("InnerProduct", &BSpline::BSplineGridRepresentation::InnerProduct)
        .def("GetGlobalGrid", &BSpline::BSplineGridRepresentation::GetGlobalGrid)
        .def("GetGlobalWeights", &BSpline::BSplineGridRepresentation::GetGlobalWeights)
        .def("SetupRepresentation", &BSpline::BSplineGridRepresentation::SetupRepresentation)
        .def("GetBSplineObject", &BSpline::BSplineGridRepresentation::GetBSplineObject)
        .def("ApplyConfigSection", &BSpline::BSplineGridRepresentation::ApplyConfigSection)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< BSpline::BSplineGridRepresentation > >();
    delete BSpline_BSplineGridRepresentation_scope;

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

    scope* BSpline_BSpline_scope = new scope(
    class_< BSpline::BSpline, boost::noncopyable >("BSpline", init<  >())
        .def_readwrite("NumberOfBSplines", &BSpline::BSpline::NumberOfBSplines)
        .def_readwrite("MaxSplineOrder", &BSpline::BSpline::MaxSplineOrder)
        .def_readwrite("ProjectionAlgorithm", &BSpline::BSpline::ProjectionAlgorithm)
        .def_readwrite("LapackAlgorithm", &BSpline::BSpline::LapackAlgorithm)
        .def("GetBreakpointSequence", &BSpline::BSpline::GetBreakpointSequence)
        .def("GetContinuitySequence", &BSpline::BSpline::GetContinuitySequence)
        .def("GetKnotSequence", &BSpline::BSpline::GetKnotSequence)
        .def("GetKnotGridIndexMap", &BSpline::BSpline::GetKnotGridIndexMap)
        .def("GetTopKnotMap", &BSpline::BSpline::GetTopKnotMap)
        .def("GetWeights", &BSpline::BSpline::GetWeights)
        .def("GetNodes", &BSpline::BSpline::GetNodes)
        .def("GetQuadratureGrid", &BSpline::BSpline::GetQuadratureGrid)
        .def("GetQuadratureGridGlobal", &BSpline::BSpline::GetQuadratureGridGlobal)
        .def("GetQuadratureWeightsGlobal", &BSpline::BSpline::GetQuadratureWeightsGlobal)
        .def("GetOverlap", &BSpline::BSpline::GetOverlap)
        .def("ResizeBreakpointSequence", &BSpline::BSpline::ResizeBreakpointSequence)
        .def("ResizeContinuitySequence", &BSpline::BSpline::ResizeContinuitySequence)
        .def("ResizeKnotSequence", &BSpline::BSpline::ResizeKnotSequence)
        .def("ResizeKnotGridIndexMap", &BSpline::BSpline::ResizeKnotGridIndexMap)
        .def("ResizeTopKnotMap", &BSpline::BSpline::ResizeTopKnotMap)
        .def("ResizeWeights", &BSpline::BSpline::ResizeWeights)
        .def("ResizeNodes", &BSpline::BSpline::ResizeNodes)
        .def("ResizeQuadratureGridGlobal", &BSpline::BSpline::ResizeQuadratureGridGlobal)
        .def("ResizeQuadratureWeightsGlobal", &BSpline::BSpline::ResizeQuadratureWeightsGlobal)
        .def("EvaluateBSpline", &BSpline::BSpline::EvaluateBSpline)
        .def("EvaluateBSplineDerivative1", &BSpline::BSpline::EvaluateBSplineDerivative1)
        .def("EvaluateBSplineDerivative2", &BSpline::BSpline::EvaluateBSplineDerivative2)
        .def("EvaluateBSplineOnGrid", &BSpline::BSpline::EvaluateBSplineOnGrid)
        .def("GetBSpline", &BSpline::BSpline::GetBSpline)
        .def("GetBSplineDerivative1", &BSpline::BSpline::GetBSplineDerivative1)
        .def("GetBSplineDerivative2", &BSpline::BSpline::GetBSplineDerivative2)
        .def("BSplineOverlapIntegral", BSpline__BSpline__BSplineOverlapIntegralint_int)
        .def("BSplineOverlapIntegral", BSpline__BSpline__BSplineOverlapIntegralblitz__Array_double_1__int_int)
        .def("ProjectOnBSpline", &BSpline::BSpline::ProjectOnBSpline)
        .def("BSplineDerivative2OverlapIntegral", &BSpline::BSpline::BSplineDerivative2OverlapIntegral)
        .def("SolveForOverlapMatrix", &BSpline::BSpline::SolveForOverlapMatrix)
        .def("BSplineGlobalOverlapIntegral", &BSpline::BSpline::BSplineGlobalOverlapIntegral_double)
        .def("BSplineGlobalOverlapIntegral", &BSpline::BSpline::BSplineGlobalOverlapIntegral_cplx)
        .def("CreateBSplineTable", &BSpline::BSpline::CreateBSplineTable)
        .def("CreateBSplineDerivativeTable", &BSpline::BSpline::CreateBSplineDerivativeTable)
        .def("SetupOverlap", &BSpline::BSpline::SetupOverlap)
        .def("ExpandFunctionInBSplines", BSpline__BSpline__ExpandFunctionInBSplinesboost__python__object)
        .def("ExpandFunctionInBSplines", BSpline__BSpline__ExpandFunctionInBSplinesblitz__Array_std__complex_double__1__blitz__Array_std__complex_double__1_)
        .def("ConstructFunctionFromBSplineExpansion", BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1_)
        .def("ConstructFunctionFromBSplineExpansion", BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1__blitz__Array_double_1_)
        .def("ConstructFunctionFromBSplineExpansion", BSpline__BSpline__ConstructFunctionFromBSplineExpansionblitz__Array_std__complex_double__1__blitz__Array_double_1__blitz__Array_std__complex_double__1_)
        .def("ComputeStartIndex", &BSpline::BSpline::ComputeStartIndex)
        .def("ComputeStopIndex", &BSpline::BSpline::ComputeStopIndex)
        .def("ComputeIndexShift", &BSpline::BSpline::ComputeIndexShift)
        .def("GetGridIndex", &BSpline::BSpline::GetGridIndex)
        .def("ScaleAndTranslate", &BSpline::BSpline::ScaleAndTranslate)
    );
    register_ptr_to_python< boost::shared_ptr< BSpline::BSpline > >();
    delete BSpline_BSpline_scope;

def("RepresentPotentialInBasisBSpline", RepresentPotentialInBasisBSpline<cplx, 1>);
def("RepresentPotentialInBasisBSpline", RepresentPotentialInBasisBSpline<cplx, 2>);
def("RepresentPotentialInBasisBSpline", RepresentPotentialInBasisBSpline<cplx, 3>);
}


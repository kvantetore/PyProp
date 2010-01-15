
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential/densematrixpotentialevaluator.h>
#include <potential/sparsematrixpotentialevaluator.h>
#include <representation/vectorrepresentation.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

cplx (VectorRepresentation::*VectorRepresentation__InnerProductconstWavefunction_1___constWavefunction_1__)(const Wavefunction<1>&, const Wavefunction<1>&)  = &VectorRepresentation::InnerProduct;

std::complex<double> (OrthogonalRepresentation::*OrthogonalRepresentation__InnerProductconstWavefunction_1___constWavefunction_1__)(const Wavefunction<1>&, const Wavefunction<1>&)  = &OrthogonalRepresentation::InnerProduct;

void (DenseMatrixPotentialEvaluator<1>::*DenseMatrixPotentialEvaluator_1___MultiplyPotentialWavefunction_1___Wavefunction_1___double)(Wavefunction<1>&, Wavefunction<1>&, double)  = &DenseMatrixPotentialEvaluator<1>::MultiplyPotential;

void (DenseMatrixPotentialEvaluator<1>::*DenseMatrixPotentialEvaluator_1___MultiplyPotentialblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1___double)(blitz::Array<std::complex<double>,1>&, blitz::Array<std::complex<double>,1>&, double)  = &DenseMatrixPotentialEvaluator<1>::MultiplyPotential;

void (DenseMatrixPotentialEvaluator<2>::*DenseMatrixPotentialEvaluator_2___MultiplyPotentialWavefunction_2___Wavefunction_2___double)(Wavefunction<2>&, Wavefunction<2>&, double)  = &DenseMatrixPotentialEvaluator<2>::MultiplyPotential;

void (DenseMatrixPotentialEvaluator<2>::*DenseMatrixPotentialEvaluator_2___MultiplyPotentialblitz__Array_std__complex_double__2___blitz__Array_std__complex_double__2___double)(blitz::Array<std::complex<double>,2>&, blitz::Array<std::complex<double>,2>&, double)  = &DenseMatrixPotentialEvaluator<2>::MultiplyPotential;


}// namespace 


// Module ======================================================================
void Export_python_vectorrepresentation()
{
    scope* VectorRepresentation_scope = new scope(
    class_< VectorRepresentation, bases< Representation<1> >  >("VectorRepresentation", init<  >())
        .def(init< const VectorRepresentation& >())
        .def_readwrite("VectorSize", &VectorRepresentation::VectorSize)
        .def_readwrite("IndexGrid", &VectorRepresentation::IndexGrid)
        .def_readwrite("Weights", &VectorRepresentation::Weights)
        .def("Copy", &VectorRepresentation::Copy)
        .def("GetGlobalGrid", &VectorRepresentation::GetGlobalGrid)
        .def("GetGlobalWeights", &VectorRepresentation::GetGlobalWeights)
        .def("GetFullShape", &VectorRepresentation::GetFullShape)
        .def("InnerProduct", VectorRepresentation__InnerProductconstWavefunction_1___constWavefunction_1__)
        .def("ApplyConfigSection", &VectorRepresentation::ApplyConfigSection)
        .def("InnerProduct", OrthogonalRepresentation__InnerProductconstWavefunction_1___constWavefunction_1__)
        .def("GetGlobalOverlapMatrix", &OrthogonalRepresentation::GetGlobalOverlapMatrix)
    );
    register_ptr_to_python< boost::shared_ptr< VectorRepresentation > >();
    delete VectorRepresentation_scope;

    class_< SparseMatrixPotentialEvaluator, boost::noncopyable >("SparseMatrixPotentialEvaluator", no_init)
        .def("SetMatrixData", &SparseMatrixPotentialEvaluator::SetMatrixData)
        .def("MultiplyPotential", &SparseMatrixPotentialEvaluator::MultiplyPotential)
        .def("CalculateExpectationValue", &SparseMatrixPotentialEvaluator::CalculateExpectationValue)
    ;

    class_< DenseMatrixPotentialEvaluator<1>, boost::noncopyable >("DenseMatrixPotentialEvaluator_1", init<  >())
        .def("GetEigenvectors", &DenseMatrixPotentialEvaluator<1>::GetEigenvectors)
        .def("GetEigenvalues", &DenseMatrixPotentialEvaluator<1>::GetEigenvalues)
        .def("GetMatrixData", &DenseMatrixPotentialEvaluator<1>::GetMatrixData)
        .def("GetAlgorithm", &DenseMatrixPotentialEvaluator<1>::GetAlgorithm)
        .def("SetAlgorithm", &DenseMatrixPotentialEvaluator<1>::SetAlgorithm)
        .def("SetMatrixData", &DenseMatrixPotentialEvaluator<1>::SetMatrixData)
        .def("DiagonalizeMatrix", &DenseMatrixPotentialEvaluator<1>::DiagonalizeMatrix)
        .def("MultiplyPotential", DenseMatrixPotentialEvaluator_1___MultiplyPotentialWavefunction_1___Wavefunction_1___double)
        .def("MultiplyPotential", DenseMatrixPotentialEvaluator_1___MultiplyPotentialblitz__Array_std__complex_double__1___blitz__Array_std__complex_double__1___double)
        .def("ApplyPotential", &DenseMatrixPotentialEvaluator<1>::ApplyPotential)
    ;

    class_< DenseMatrixPotentialEvaluator<2>, boost::noncopyable >("DenseMatrixPotentialEvaluator_2", init<  >())
        .def("GetEigenvectors", &DenseMatrixPotentialEvaluator<2>::GetEigenvectors)
        .def("GetEigenvalues", &DenseMatrixPotentialEvaluator<2>::GetEigenvalues)
        .def("GetMatrixData", &DenseMatrixPotentialEvaluator<2>::GetMatrixData)
        .def("GetAlgorithm", &DenseMatrixPotentialEvaluator<2>::GetAlgorithm)
        .def("SetAlgorithm", &DenseMatrixPotentialEvaluator<2>::SetAlgorithm)
        .def("SetMatrixData", &DenseMatrixPotentialEvaluator<2>::SetMatrixData)
        .def("DiagonalizeMatrix", &DenseMatrixPotentialEvaluator<2>::DiagonalizeMatrix)
        .def("MultiplyPotential", DenseMatrixPotentialEvaluator_2___MultiplyPotentialWavefunction_2___Wavefunction_2___double)
        .def("MultiplyPotential", DenseMatrixPotentialEvaluator_2___MultiplyPotentialblitz__Array_std__complex_double__2___blitz__Array_std__complex_double__2___double)
        .def("ApplyPotential", &DenseMatrixPotentialEvaluator<2>::ApplyPotential)
    ;

    class_< DenseMatrixPotentialEvaluatorOld, boost::noncopyable >("DenseMatrixPotentialEvaluatorOld", no_init)
        .def("SetMatrixData", &DenseMatrixPotentialEvaluatorOld::SetMatrixData)
        .def("MultiplyPotential", &DenseMatrixPotentialEvaluatorOld::MultiplyPotential)
        .def("CalculateExpectationValue", &DenseMatrixPotentialEvaluatorOld::CalculateExpectationValue)
    ;

}


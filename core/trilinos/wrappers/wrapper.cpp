
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <../pyprop_epetra.cpp>
#include <epetrapotential.h>
#include <ifpackwrapper.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libtrilinos)
{
    class_< IfpackRadialPreconditioner<1> >("IfpackRadialPreconditioner_1", init<  >())
        .def(init< const IfpackRadialPreconditioner<1>& >())
        .def("ApplyConfigSection", &IfpackRadialPreconditioner<1>::ApplyConfigSection)
        .def("Setup", &IfpackRadialPreconditioner<1>::Setup)
        .def("Solve", &IfpackRadialPreconditioner<1>::Solve)
    ;

    class_< IfpackRadialPreconditioner<2> >("IfpackRadialPreconditioner_2", init<  >())
        .def(init< const IfpackRadialPreconditioner<2>& >())
        .def("ApplyConfigSection", &IfpackRadialPreconditioner<2>::ApplyConfigSection)
        .def("Setup", &IfpackRadialPreconditioner<2>::Setup)
        .def("Solve", &IfpackRadialPreconditioner<2>::Solve)
    ;

    class_< IfpackRadialPreconditioner<3> >("IfpackRadialPreconditioner_3", init<  >())
        .def(init< const IfpackRadialPreconditioner<3>& >())
        .def("ApplyConfigSection", &IfpackRadialPreconditioner<3>::ApplyConfigSection)
        .def("Setup", &IfpackRadialPreconditioner<3>::Setup)
        .def("Solve", &IfpackRadialPreconditioner<3>::Solve)
    ;

    class_< EpetraPotential<1> >("EpetraPotential_1", init<  >())
        .def(init< const EpetraPotential<1>& >())
        .def("Setup", &EpetraPotential<1>::Setup)
        .def("AddTensorPotentialData", &EpetraPotential<1>::AddTensorPotentialData)
        .def("GlobalAssemble", &EpetraPotential<1>::GlobalAssemble)
        .def("Filled", &EpetraPotential<1>::Filled)
        .def("Multiply", &EpetraPotential<1>::Multiply)
        .def("MultiplyPotential", &EpetraPotential<1>::MultiplyPotential)
        .def("NumMyNonzeros", &EpetraPotential<1>::NumMyNonzeros)
        .def("NumGlobalRows", &EpetraPotential<1>::NumGlobalRows)
        .def("NumGlobalCols", &EpetraPotential<1>::NumGlobalCols)
        .def("StorageOptimizied", &EpetraPotential<1>::StorageOptimizied)
    ;

    class_< EpetraPotential<2> >("EpetraPotential_2", init<  >())
        .def(init< const EpetraPotential<2>& >())
        .def("Setup", &EpetraPotential<2>::Setup)
        .def("AddTensorPotentialData", &EpetraPotential<2>::AddTensorPotentialData)
        .def("GlobalAssemble", &EpetraPotential<2>::GlobalAssemble)
        .def("Filled", &EpetraPotential<2>::Filled)
        .def("Multiply", &EpetraPotential<2>::Multiply)
        .def("MultiplyPotential", &EpetraPotential<2>::MultiplyPotential)
        .def("NumMyNonzeros", &EpetraPotential<2>::NumMyNonzeros)
        .def("NumGlobalRows", &EpetraPotential<2>::NumGlobalRows)
        .def("NumGlobalCols", &EpetraPotential<2>::NumGlobalCols)
        .def("StorageOptimizied", &EpetraPotential<2>::StorageOptimizied)
    ;

    class_< EpetraPotential<3> >("EpetraPotential_3", init<  >())
        .def(init< const EpetraPotential<3>& >())
        .def("Setup", &EpetraPotential<3>::Setup)
        .def("AddTensorPotentialData", &EpetraPotential<3>::AddTensorPotentialData)
        .def("GlobalAssemble", &EpetraPotential<3>::GlobalAssemble)
        .def("Filled", &EpetraPotential<3>::Filled)
        .def("Multiply", &EpetraPotential<3>::Multiply)
        .def("MultiplyPotential", &EpetraPotential<3>::MultiplyPotential)
        .def("NumMyNonzeros", &EpetraPotential<3>::NumMyNonzeros)
        .def("NumGlobalRows", &EpetraPotential<3>::NumGlobalRows)
        .def("NumGlobalCols", &EpetraPotential<3>::NumGlobalCols)
        .def("StorageOptimizied", &EpetraPotential<3>::StorageOptimizied)
    ;

    class_< EpetraPotential<4> >("EpetraPotential_4", init<  >())
        .def(init< const EpetraPotential<4>& >())
        .def("Setup", &EpetraPotential<4>::Setup)
        .def("AddTensorPotentialData", &EpetraPotential<4>::AddTensorPotentialData)
        .def("GlobalAssemble", &EpetraPotential<4>::GlobalAssemble)
        .def("Filled", &EpetraPotential<4>::Filled)
        .def("Multiply", &EpetraPotential<4>::Multiply)
        .def("MultiplyPotential", &EpetraPotential<4>::MultiplyPotential)
        .def("NumMyNonzeros", &EpetraPotential<4>::NumMyNonzeros)
        .def("NumGlobalRows", &EpetraPotential<4>::NumGlobalRows)
        .def("NumGlobalCols", &EpetraPotential<4>::NumGlobalCols)
        .def("StorageOptimizied", &EpetraPotential<4>::StorageOptimizied)
    ;

def("CreateTensorPotentialEpetraMatrix_1", CreateTensorPotentialEpetraMatrix<1>);
def("CreateTensorPotentialEpetraMatrix_2", CreateTensorPotentialEpetraMatrix<2>);
def("CreateTensorPotentialEpetraMatrix_3", CreateTensorPotentialEpetraMatrix<3>);
def("CreateTensorPotentialEpetraMatrix_4", CreateTensorPotentialEpetraMatrix<4>);
}


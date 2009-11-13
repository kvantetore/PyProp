
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
    class_< IfpackRadialPreconditioner<1>, boost::noncopyable >("IfpackRadialPreconditioner_1", init<  >())
        .def("ApplyConfigSection", &IfpackRadialPreconditioner<1>::ApplyConfigSection)
        .def("Setup", &IfpackRadialPreconditioner<1>::Setup)
        .def("Solve", &IfpackRadialPreconditioner<1>::Solve)
    ;

    class_< IfpackRadialPreconditioner<2>, boost::noncopyable >("IfpackRadialPreconditioner_2", init<  >())
        .def("ApplyConfigSection", &IfpackRadialPreconditioner<2>::ApplyConfigSection)
        .def("Setup", &IfpackRadialPreconditioner<2>::Setup)
        .def("Solve", &IfpackRadialPreconditioner<2>::Solve)
    ;

    class_< IfpackRadialPreconditioner<3>, boost::noncopyable >("IfpackRadialPreconditioner_3", init<  >())
        .def("ApplyConfigSection", &IfpackRadialPreconditioner<3>::ApplyConfigSection)
        .def("Setup", &IfpackRadialPreconditioner<3>::Setup)
        .def("Solve", &IfpackRadialPreconditioner<3>::Solve)
    ;

    class_< EpetraPotential<1>, boost::noncopyable >("EpetraPotential_1", init<  >())
        .def("Setup", &EpetraPotential<1>::Setup)
        .def("AddTensorPotentialData", &EpetraPotential<1>::AddTensorPotentialData)
        .def("GlobalAssemble", &EpetraPotential<1>::GlobalAssemble)
        .def("Filled", &EpetraPotential<1>::Filled)
        .def("Multiply", &EpetraPotential<1>::Multiply)
    ;

    class_< EpetraPotential<2>, boost::noncopyable >("EpetraPotential_2", init<  >())
        .def("Setup", &EpetraPotential<2>::Setup)
        .def("AddTensorPotentialData", &EpetraPotential<2>::AddTensorPotentialData)
        .def("GlobalAssemble", &EpetraPotential<2>::GlobalAssemble)
        .def("Filled", &EpetraPotential<2>::Filled)
        .def("Multiply", &EpetraPotential<2>::Multiply)
    ;

    class_< EpetraPotential<3>, boost::noncopyable >("EpetraPotential_3", init<  >())
        .def("Setup", &EpetraPotential<3>::Setup)
        .def("AddTensorPotentialData", &EpetraPotential<3>::AddTensorPotentialData)
        .def("GlobalAssemble", &EpetraPotential<3>::GlobalAssemble)
        .def("Filled", &EpetraPotential<3>::Filled)
        .def("Multiply", &EpetraPotential<3>::Multiply)
    ;

    class_< EpetraPotential<4>, boost::noncopyable >("EpetraPotential_4", init<  >())
        .def("Setup", &EpetraPotential<4>::Setup)
        .def("AddTensorPotentialData", &EpetraPotential<4>::AddTensorPotentialData)
        .def("GlobalAssemble", &EpetraPotential<4>::GlobalAssemble)
        .def("Filled", &EpetraPotential<4>::Filled)
        .def("Multiply", &EpetraPotential<4>::Multiply)
    ;

def("CreateTensorPotentialEpetraMatrix_1", CreateTensorPotentialEpetraMatrix<1>);
def("CreateTensorPotentialEpetraMatrix_2", CreateTensorPotentialEpetraMatrix<2>);
def("CreateTensorPotentialEpetraMatrix_3", CreateTensorPotentialEpetraMatrix<3>);
def("CreateTensorPotentialEpetraMatrix_4", CreateTensorPotentialEpetraMatrix<4>);
}


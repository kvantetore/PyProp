
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <../pyprop_epetra.cpp>
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

def("CreateTensorPotentialEpetraMatrix_1", CreateTensorPotentialEpetraMatrix<1>);
def("CreateTensorPotentialEpetraMatrix_2", CreateTensorPotentialEpetraMatrix<2>);
def("CreateTensorPotentialEpetraMatrix_3", CreateTensorPotentialEpetraMatrix<3>);
def("CreateTensorPotentialEpetraMatrix_4", CreateTensorPotentialEpetraMatrix<4>);
}


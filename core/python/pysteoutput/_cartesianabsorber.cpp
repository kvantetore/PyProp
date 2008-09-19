
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential/cartesianabsorber.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_cartesianabsorber()
{
    class_< CartesianAbsorbingPotential<1> >("CartesianAbsorbingPotential_1", init<  >())
        .def(init< const CartesianAbsorbingPotential<1>& >())
        .def("ApplyPotential", &CartesianAbsorbingPotential<1>::ApplyPotential)
    ;

    class_< CartesianAbsorbingPotential<2> >("CartesianAbsorbingPotential_2", init<  >())
        .def(init< const CartesianAbsorbingPotential<2>& >())
        .def("ApplyPotential", &CartesianAbsorbingPotential<2>::ApplyPotential)
    ;

    class_< CartesianAbsorbingPotential<3> >("CartesianAbsorbingPotential_3", init<  >())
        .def(init< const CartesianAbsorbingPotential<3>& >())
        .def("ApplyPotential", &CartesianAbsorbingPotential<3>::ApplyPotential)
    ;

    class_< CartesianAbsorbingPotential<4> >("CartesianAbsorbingPotential_4", init<  >())
        .def(init< const CartesianAbsorbingPotential<4>& >())
        .def("ApplyPotential", &CartesianAbsorbingPotential<4>::ApplyPotential)
    ;

    class_< CartesianBoundaryAbsorbingPotential<1> >("CartesianBoundaryAbsorbingPotential_1", init<  >())
        .def(init< const CartesianBoundaryAbsorbingPotential<1>& >())
        .def_readwrite("AbsorbtionWidth", &CartesianBoundaryAbsorbingPotential<1>::AbsorbtionWidth)
        .def("ApplyConfigSection", &CartesianBoundaryAbsorbingPotential<1>::ApplyConfigSection)
        .def("ApplyPotential", &CartesianBoundaryAbsorbingPotential<1>::ApplyPotential)
    ;

    class_< CartesianBoundaryAbsorbingPotential<2> >("CartesianBoundaryAbsorbingPotential_2", init<  >())
        .def(init< const CartesianBoundaryAbsorbingPotential<2>& >())
        .def_readwrite("AbsorbtionWidth", &CartesianBoundaryAbsorbingPotential<2>::AbsorbtionWidth)
        .def("ApplyConfigSection", &CartesianBoundaryAbsorbingPotential<2>::ApplyConfigSection)
        .def("ApplyPotential", &CartesianBoundaryAbsorbingPotential<2>::ApplyPotential)
    ;

    class_< CartesianBoundaryAbsorbingPotential<3> >("CartesianBoundaryAbsorbingPotential_3", init<  >())
        .def(init< const CartesianBoundaryAbsorbingPotential<3>& >())
        .def_readwrite("AbsorbtionWidth", &CartesianBoundaryAbsorbingPotential<3>::AbsorbtionWidth)
        .def("ApplyConfigSection", &CartesianBoundaryAbsorbingPotential<3>::ApplyConfigSection)
        .def("ApplyPotential", &CartesianBoundaryAbsorbingPotential<3>::ApplyPotential)
    ;

    class_< CartesianBoundaryAbsorbingPotential<4> >("CartesianBoundaryAbsorbingPotential_4", init<  >())
        .def(init< const CartesianBoundaryAbsorbingPotential<4>& >())
        .def_readwrite("AbsorbtionWidth", &CartesianBoundaryAbsorbingPotential<4>::AbsorbtionWidth)
        .def("ApplyConfigSection", &CartesianBoundaryAbsorbingPotential<4>::ApplyConfigSection)
        .def("ApplyPotential", &CartesianBoundaryAbsorbingPotential<4>::ApplyPotential)
    ;

}


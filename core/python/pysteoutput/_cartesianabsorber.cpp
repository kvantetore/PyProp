
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
    class_< CartesianAbsorbingPotential<1>, boost::noncopyable >("CartesianAbsorbingPotential_1", no_init)
        .def("ApplyPotential", &CartesianAbsorbingPotential<1>::ApplyPotential)
    ;

    class_< CartesianAbsorbingPotential<2>, boost::noncopyable >("CartesianAbsorbingPotential_2", no_init)
        .def("ApplyPotential", &CartesianAbsorbingPotential<2>::ApplyPotential)
    ;

    class_< CartesianAbsorbingPotential<3>, boost::noncopyable >("CartesianAbsorbingPotential_3", no_init)
        .def("ApplyPotential", &CartesianAbsorbingPotential<3>::ApplyPotential)
    ;

    class_< CartesianAbsorbingPotential<4>, boost::noncopyable >("CartesianAbsorbingPotential_4", no_init)
        .def("ApplyPotential", &CartesianAbsorbingPotential<4>::ApplyPotential)
    ;

    class_< CartesianBoundaryAbsorbingPotential<1>, boost::noncopyable >("CartesianBoundaryAbsorbingPotential_1", no_init)
        .def_readwrite("AbsorbtionWidth", &CartesianBoundaryAbsorbingPotential<1>::AbsorbtionWidth)
        .def("ApplyConfigSection", &CartesianBoundaryAbsorbingPotential<1>::ApplyConfigSection)
        .def("ApplyPotential", &CartesianBoundaryAbsorbingPotential<1>::ApplyPotential)
    ;

    class_< CartesianBoundaryAbsorbingPotential<2>, boost::noncopyable >("CartesianBoundaryAbsorbingPotential_2", no_init)
        .def_readwrite("AbsorbtionWidth", &CartesianBoundaryAbsorbingPotential<2>::AbsorbtionWidth)
        .def("ApplyConfigSection", &CartesianBoundaryAbsorbingPotential<2>::ApplyConfigSection)
        .def("ApplyPotential", &CartesianBoundaryAbsorbingPotential<2>::ApplyPotential)
    ;

    class_< CartesianBoundaryAbsorbingPotential<3>, boost::noncopyable >("CartesianBoundaryAbsorbingPotential_3", no_init)
        .def_readwrite("AbsorbtionWidth", &CartesianBoundaryAbsorbingPotential<3>::AbsorbtionWidth)
        .def("ApplyConfigSection", &CartesianBoundaryAbsorbingPotential<3>::ApplyConfigSection)
        .def("ApplyPotential", &CartesianBoundaryAbsorbingPotential<3>::ApplyPotential)
    ;

    class_< CartesianBoundaryAbsorbingPotential<4>, boost::noncopyable >("CartesianBoundaryAbsorbingPotential_4", no_init)
        .def_readwrite("AbsorbtionWidth", &CartesianBoundaryAbsorbingPotential<4>::AbsorbtionWidth)
        .def("ApplyConfigSection", &CartesianBoundaryAbsorbingPotential<4>::ApplyConfigSection)
        .def("ApplyPotential", &CartesianBoundaryAbsorbingPotential<4>::ApplyPotential)
    ;

}



// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential/sphericalabsorber.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_sphericalabsorber()
{
    class_< SphericalAbsorbingPotential<2>, boost::noncopyable >("SphericalAbsorbingPotential_2", no_init)
        .def("ApplyConfigSection", &SphericalAbsorbingPotential<2>::ApplyConfigSection)
        .def("ApplyPotential", &SphericalAbsorbingPotential<2>::ApplyPotential)
        .def("MultiplyPotential", &SphericalAbsorbingPotential<2>::MultiplyPotential)
    ;

    class_< SphericalAbsorbingPotential<3>, boost::noncopyable >("SphericalAbsorbingPotential_3", no_init)
        .def("ApplyConfigSection", &SphericalAbsorbingPotential<3>::ApplyConfigSection)
        .def("ApplyPotential", &SphericalAbsorbingPotential<3>::ApplyPotential)
        .def("MultiplyPotential", &SphericalAbsorbingPotential<3>::MultiplyPotential)
    ;

}



// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <representation/spherical/lmrange.h>
#include <representation/spherical/omegarange.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_sphericalrange()
{
    scope* OmegaRange_scope = new scope(
    class_< OmegaRange, boost::noncopyable >("OmegaRange", init<  >())
        .def_readwrite("Type", &OmegaRange::Type)
        .def_readwrite("MaxL", &OmegaRange::MaxL)
        .def("SetupRange", &OmegaRange::SetupRange)
        .def("Count", &OmegaRange::Count)
        .def("GetIndexGrid", &OmegaRange::GetIndexGrid, return_value_policy< return_by_value >())
        .def("GetOmegaGrid", &OmegaRange::GetOmegaGrid, return_value_policy< return_by_value >())
        .def("GetWeights", &OmegaRange::GetWeights, return_value_policy< copy_const_reference >())
    );

    enum_< OmegaRange::GridType >("GridType")
        .value("SloanWomersley", OmegaRange::SloanWomersley)
        .value("Equidistant", OmegaRange::Equidistant)
        .value("Gauss", OmegaRange::Gauss)
    ;

    delete OmegaRange_scope;

    class_< LmRange, boost::noncopyable >("LmRange", init<  >())
        .def(init< int >())
        .def_readwrite("MaxL", &LmRange::MaxL)
        .def("Count", &LmRange::Count)
        .def("GetIndexGrid", &LmRange::GetIndexGrid, return_value_policy< return_by_value >())
        .def("GetLmGrid", &LmRange::GetLmGrid, return_value_policy< copy_const_reference >())
        .def("GetWeights", &LmRange::GetWeights)
        .def("MapLmIndex", &LmRange::MapLmIndex)
        .def("GetL", &LmRange::GetL)
        .def("GetM", &LmRange::GetM)
        .def("MapLmIndexInverse", &LmRange::MapLmIndexInverse)
        .staticmethod("MapLmIndexInverse")
        .staticmethod("GetL")
        .staticmethod("GetM")
        .staticmethod("MapLmIndex")
    ;

}


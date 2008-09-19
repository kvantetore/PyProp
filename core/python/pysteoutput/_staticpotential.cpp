
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential/staticpotential.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_staticpotential()
{
    scope* StaticPotential_1_scope = new scope(
    class_< StaticPotential<1> >("StaticPotential_1", init<  >())
        .def(init< const StaticPotential<1>& >())
        .def("InitializePotential", &StaticPotential<1>::InitializePotential)
        .def("GetPotentialData", &StaticPotential<1>::GetPotentialData)
        .def("GetStorageModel", &StaticPotential<1>::GetStorageModel)
        .def("ApplyPotential", &StaticPotential<1>::ApplyPotential)
        .def("MultiplyPotential", &StaticPotential<1>::MultiplyPotential)
    );

    enum_< StaticPotential<1>::StorageModel >("StorageModel")
        .value("StorageExpValue", StaticPotential<1>::StorageExpValue)
        .value("StorageValue", StaticPotential<1>::StorageValue)
    ;

    delete StaticPotential_1_scope;

    scope* StaticPotential_2_scope = new scope(
    class_< StaticPotential<2> >("StaticPotential_2", init<  >())
        .def(init< const StaticPotential<2>& >())
        .def("InitializePotential", &StaticPotential<2>::InitializePotential)
        .def("GetPotentialData", &StaticPotential<2>::GetPotentialData)
        .def("GetStorageModel", &StaticPotential<2>::GetStorageModel)
        .def("ApplyPotential", &StaticPotential<2>::ApplyPotential)
        .def("MultiplyPotential", &StaticPotential<2>::MultiplyPotential)
    );

    enum_< StaticPotential<2>::StorageModel >("StorageModel")
        .value("StorageExpValue", StaticPotential<2>::StorageExpValue)
        .value("StorageValue", StaticPotential<2>::StorageValue)
    ;

    delete StaticPotential_2_scope;

    scope* StaticPotential_3_scope = new scope(
    class_< StaticPotential<3> >("StaticPotential_3", init<  >())
        .def(init< const StaticPotential<3>& >())
        .def("InitializePotential", &StaticPotential<3>::InitializePotential)
        .def("GetPotentialData", &StaticPotential<3>::GetPotentialData)
        .def("GetStorageModel", &StaticPotential<3>::GetStorageModel)
        .def("ApplyPotential", &StaticPotential<3>::ApplyPotential)
        .def("MultiplyPotential", &StaticPotential<3>::MultiplyPotential)
    );

    enum_< StaticPotential<3>::StorageModel >("StorageModel")
        .value("StorageExpValue", StaticPotential<3>::StorageExpValue)
        .value("StorageValue", StaticPotential<3>::StorageValue)
    ;

    delete StaticPotential_3_scope;

    scope* StaticPotential_4_scope = new scope(
    class_< StaticPotential<4> >("StaticPotential_4", init<  >())
        .def(init< const StaticPotential<4>& >())
        .def("InitializePotential", &StaticPotential<4>::InitializePotential)
        .def("GetPotentialData", &StaticPotential<4>::GetPotentialData)
        .def("GetStorageModel", &StaticPotential<4>::GetStorageModel)
        .def("ApplyPotential", &StaticPotential<4>::ApplyPotential)
        .def("MultiplyPotential", &StaticPotential<4>::MultiplyPotential)
    );

    enum_< StaticPotential<4>::StorageModel >("StorageModel")
        .value("StorageExpValue", StaticPotential<4>::StorageExpValue)
        .value("StorageValue", StaticPotential<4>::StorageValue)
    ;

    delete StaticPotential_4_scope;

}



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
    class_< StaticPotential<1>, boost::noncopyable >("StaticPotential_1", init<  >())
        .def("UseStorageValue", &StaticPotential<1>::UseStorageValue)
        .def("UseStorageExpValue", &StaticPotential<1>::UseStorageExpValue)
        .def("InitializePotential", &StaticPotential<1>::InitializePotential)
        .def("GetPotentialData", &StaticPotential<1>::GetPotentialData)
        .def("GetPotentialDataExp", &StaticPotential<1>::GetPotentialDataExp)
        .def("GetStorageModel", &StaticPotential<1>::GetStorageModel)
        .def("ApplyPotential", &StaticPotential<1>::ApplyPotential)
        .def("MultiplyPotential", &StaticPotential<1>::MultiplyPotential)
    );

    enum_< StaticPotential<1>::StorageModel >("StorageModel")
        .value("StorageExpValue", StaticPotential<1>::StorageExpValue)
        .value("StorageBoth", StaticPotential<1>::StorageBoth)
        .value("StorageValue", StaticPotential<1>::StorageValue)
    ;

    delete StaticPotential_1_scope;

    scope* StaticPotential_2_scope = new scope(
    class_< StaticPotential<2>, boost::noncopyable >("StaticPotential_2", init<  >())
        .def("UseStorageValue", &StaticPotential<2>::UseStorageValue)
        .def("UseStorageExpValue", &StaticPotential<2>::UseStorageExpValue)
        .def("InitializePotential", &StaticPotential<2>::InitializePotential)
        .def("GetPotentialData", &StaticPotential<2>::GetPotentialData)
        .def("GetPotentialDataExp", &StaticPotential<2>::GetPotentialDataExp)
        .def("GetStorageModel", &StaticPotential<2>::GetStorageModel)
        .def("ApplyPotential", &StaticPotential<2>::ApplyPotential)
        .def("MultiplyPotential", &StaticPotential<2>::MultiplyPotential)
    );

    enum_< StaticPotential<2>::StorageModel >("StorageModel")
        .value("StorageExpValue", StaticPotential<2>::StorageExpValue)
        .value("StorageBoth", StaticPotential<2>::StorageBoth)
        .value("StorageValue", StaticPotential<2>::StorageValue)
    ;

    delete StaticPotential_2_scope;

    scope* StaticPotential_3_scope = new scope(
    class_< StaticPotential<3>, boost::noncopyable >("StaticPotential_3", init<  >())
        .def("UseStorageValue", &StaticPotential<3>::UseStorageValue)
        .def("UseStorageExpValue", &StaticPotential<3>::UseStorageExpValue)
        .def("InitializePotential", &StaticPotential<3>::InitializePotential)
        .def("GetPotentialData", &StaticPotential<3>::GetPotentialData)
        .def("GetPotentialDataExp", &StaticPotential<3>::GetPotentialDataExp)
        .def("GetStorageModel", &StaticPotential<3>::GetStorageModel)
        .def("ApplyPotential", &StaticPotential<3>::ApplyPotential)
        .def("MultiplyPotential", &StaticPotential<3>::MultiplyPotential)
    );

    enum_< StaticPotential<3>::StorageModel >("StorageModel")
        .value("StorageExpValue", StaticPotential<3>::StorageExpValue)
        .value("StorageBoth", StaticPotential<3>::StorageBoth)
        .value("StorageValue", StaticPotential<3>::StorageValue)
    ;

    delete StaticPotential_3_scope;

    scope* StaticPotential_4_scope = new scope(
    class_< StaticPotential<4>, boost::noncopyable >("StaticPotential_4", init<  >())
        .def("UseStorageValue", &StaticPotential<4>::UseStorageValue)
        .def("UseStorageExpValue", &StaticPotential<4>::UseStorageExpValue)
        .def("InitializePotential", &StaticPotential<4>::InitializePotential)
        .def("GetPotentialData", &StaticPotential<4>::GetPotentialData)
        .def("GetPotentialDataExp", &StaticPotential<4>::GetPotentialDataExp)
        .def("GetStorageModel", &StaticPotential<4>::GetStorageModel)
        .def("ApplyPotential", &StaticPotential<4>::ApplyPotential)
        .def("MultiplyPotential", &StaticPotential<4>::MultiplyPotential)
    );

    enum_< StaticPotential<4>::StorageModel >("StorageModel")
        .value("StorageExpValue", StaticPotential<4>::StorageExpValue)
        .value("StorageBoth", StaticPotential<4>::StorageBoth)
        .value("StorageValue", StaticPotential<4>::StorageValue)
    ;

    delete StaticPotential_4_scope;

}



// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential/combinedabsorber.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct AbsorberModel_Wrapper: AbsorberModel
{
    AbsorberModel_Wrapper(PyObject* py_self_, const AbsorberModel& p0):
        AbsorberModel(p0), py_self(py_self_) {}

    AbsorberModel_Wrapper(PyObject* py_self_):
        AbsorberModel(), py_self(py_self_) {}

    void ApplyConfigSection(const ConfigSection& p0) {
        call_method< void >(py_self, "ApplyConfigSection", p0);
    }

    void default_ApplyConfigSection(const ConfigSection& p0) {
        AbsorberModel::ApplyConfigSection(p0);
    }

    void SetupStep(blitz::Array<double,1> p0) {
        call_method< void >(py_self, "SetupStep", p0);
    }

    blitz::Array<double,1> GetScaling() {
        return call_method< blitz::Array<double,1> >(py_self, "GetScaling");
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
void Export_python_combinedabsorber()
{
    scope* AbsorberModel_scope = new scope(
    class_< AbsorberModel, boost::noncopyable, AbsorberModel_Wrapper >("AbsorberModel", init<  >())
        .def("ApplyConfigSection", &AbsorberModel::ApplyConfigSection, &AbsorberModel_Wrapper::default_ApplyConfigSection)
        .def("SetupStep", pure_virtual(&AbsorberModel::SetupStep))
        .def("GetScaling", pure_virtual(&AbsorberModel::GetScaling))
    );
    register_ptr_to_python< boost::shared_ptr< AbsorberModel > >();
    delete AbsorberModel_scope;

    class_< CombinedAbsorberPotential<1> >("CombinedAbsorberPotential_1", init<  >())
        .def(init< const CombinedAbsorberPotential<1>& >())
        .def("AddAbsorber", &CombinedAbsorberPotential<1>::AddAbsorber)
        .def("ApplyPotential", &CombinedAbsorberPotential<1>::ApplyPotential)
    ;

    class_< CombinedAbsorberPotential<2> >("CombinedAbsorberPotential_2", init<  >())
        .def(init< const CombinedAbsorberPotential<2>& >())
        .def("AddAbsorber", &CombinedAbsorberPotential<2>::AddAbsorber)
        .def("ApplyPotential", &CombinedAbsorberPotential<2>::ApplyPotential)
    ;

    class_< CombinedAbsorberPotential<3> >("CombinedAbsorberPotential_3", init<  >())
        .def(init< const CombinedAbsorberPotential<3>& >())
        .def("AddAbsorber", &CombinedAbsorberPotential<3>::AddAbsorber)
        .def("ApplyPotential", &CombinedAbsorberPotential<3>::ApplyPotential)
    ;

    class_< CombinedAbsorberPotential<4> >("CombinedAbsorberPotential_4", init<  >())
        .def(init< const CombinedAbsorberPotential<4>& >())
        .def("AddAbsorber", &CombinedAbsorberPotential<4>::AddAbsorber)
        .def("ApplyPotential", &CombinedAbsorberPotential<4>::ApplyPotential)
    ;

}


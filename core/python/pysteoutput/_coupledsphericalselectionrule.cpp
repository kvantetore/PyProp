
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <tensorpotential/coupledsphericalselectionrule.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct CoupledSphericalSelectionRule_Wrapper: CoupledSphericalSelectionRule
{
    CoupledSphericalSelectionRule_Wrapper(PyObject* py_self_, const CoupledSphericalSelectionRule& p0):
        CoupledSphericalSelectionRule(p0), py_self(py_self_) {}

    CoupledSphericalSelectionRule_Wrapper(PyObject* py_self_):
        CoupledSphericalSelectionRule(), py_self(py_self_) {}

    bool SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    PyObject* py_self;
};

struct CoupledSphericalSelectionRuleR12_Wrapper: CoupledSphericalSelectionRuleR12
{
    CoupledSphericalSelectionRuleR12_Wrapper(PyObject* py_self_):
        CoupledSphericalSelectionRuleR12(), py_self(py_self_) {}

    CoupledSphericalSelectionRuleR12_Wrapper(PyObject* py_self_, int p0):
        CoupledSphericalSelectionRuleR12(p0), py_self(py_self_) {}

    bool SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return CoupledSphericalSelectionRuleR12::SelectionRule(p0, p1);
    }

    PyObject* py_self;
};

struct CoupledSphericalSelectionRuleR12Old_Wrapper: CoupledSphericalSelectionRuleR12Old
{
    CoupledSphericalSelectionRuleR12Old_Wrapper(PyObject* py_self_):
        CoupledSphericalSelectionRuleR12Old(), py_self(py_self_) {}

    bool SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return CoupledSphericalSelectionRuleR12Old::SelectionRule(p0, p1);
    }

    PyObject* py_self;
};

struct CoupledSphericalSelectionRuleLinearPolarizedField_Wrapper: CoupledSphericalSelectionRuleLinearPolarizedField
{
    CoupledSphericalSelectionRuleLinearPolarizedField_Wrapper(PyObject* py_self_):
        CoupledSphericalSelectionRuleLinearPolarizedField(), py_self(py_self_) {}

    bool SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return CoupledSphericalSelectionRuleLinearPolarizedField::SelectionRule(p0, p1);
    }

    PyObject* py_self;
};

struct CoupledSphericalSelectionRuleDiagonal_Wrapper: CoupledSphericalSelectionRuleDiagonal
{
    CoupledSphericalSelectionRuleDiagonal_Wrapper(PyObject* py_self_):
        CoupledSphericalSelectionRuleDiagonal(), py_self(py_self_) {}

    bool SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return CoupledSphericalSelectionRuleDiagonal::SelectionRule(p0, p1);
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
void Export_python_coupledsphericalselectionrule()
{
    class_< CoupledSphericalSelectionRule, boost::noncopyable, CoupledSphericalSelectionRule_Wrapper >("CoupledSphericalSelectionRule", init<  >())
        .def("SelectionRule", pure_virtual(&CoupledSphericalSelectionRule::SelectionRule))
        .def("GetBasisPairs", &CoupledSphericalSelectionRule::GetBasisPairs)
    ;

    class_< CoupledSphericalSelectionRuleR12, bases< CoupledSphericalSelectionRule > , boost::noncopyable, CoupledSphericalSelectionRuleR12_Wrapper >("CoupledSphericalSelectionRuleR12", init<  >())
        .def(init< int >())
        .def_readwrite("cg", &CoupledSphericalSelectionRuleR12::cg)
        .def("SelectionRule", (bool (CoupledSphericalSelectionRuleR12::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&) )&CoupledSphericalSelectionRuleR12::SelectionRule, (bool (CoupledSphericalSelectionRuleR12_Wrapper::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&))&CoupledSphericalSelectionRuleR12_Wrapper::default_SelectionRule)
    ;

    class_< CoupledSphericalSelectionRuleR12Old, bases< CoupledSphericalSelectionRule > , boost::noncopyable, CoupledSphericalSelectionRuleR12Old_Wrapper >("CoupledSphericalSelectionRuleR12Old", init<  >())
        .def_readwrite("cg", &CoupledSphericalSelectionRuleR12Old::cg)
        .def("SelectionRule", (bool (CoupledSphericalSelectionRuleR12Old::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&) )&CoupledSphericalSelectionRuleR12Old::SelectionRule, (bool (CoupledSphericalSelectionRuleR12Old_Wrapper::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&))&CoupledSphericalSelectionRuleR12Old_Wrapper::default_SelectionRule)
    ;

    class_< CoupledSphericalSelectionRuleLinearPolarizedField, bases< CoupledSphericalSelectionRule > , boost::noncopyable, CoupledSphericalSelectionRuleLinearPolarizedField_Wrapper >("CoupledSphericalSelectionRuleLinearPolarizedField", init<  >())
        .def("SelectionRule", (bool (CoupledSphericalSelectionRuleLinearPolarizedField::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&) )&CoupledSphericalSelectionRuleLinearPolarizedField::SelectionRule, (bool (CoupledSphericalSelectionRuleLinearPolarizedField_Wrapper::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&))&CoupledSphericalSelectionRuleLinearPolarizedField_Wrapper::default_SelectionRule)
    ;

    class_< CoupledSphericalSelectionRuleDiagonal, bases< CoupledSphericalSelectionRule > , boost::noncopyable, CoupledSphericalSelectionRuleDiagonal_Wrapper >("CoupledSphericalSelectionRuleDiagonal", init<  >())
        .def("SelectionRule", (bool (CoupledSphericalSelectionRuleDiagonal::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&) )&CoupledSphericalSelectionRuleDiagonal::SelectionRule, (bool (CoupledSphericalSelectionRuleDiagonal_Wrapper::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&))&CoupledSphericalSelectionRuleDiagonal_Wrapper::default_SelectionRule)
    ;

}


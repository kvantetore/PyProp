
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
    CoupledSphericalSelectionRuleR12_Wrapper(PyObject* py_self_, const CoupledSphericalSelectionRuleR12& p0):
        CoupledSphericalSelectionRuleR12(p0), py_self(py_self_) {}

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
    CoupledSphericalSelectionRuleR12Old_Wrapper(PyObject* py_self_, const CoupledSphericalSelectionRuleR12Old& p0):
        CoupledSphericalSelectionRuleR12Old(p0), py_self(py_self_) {}

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
    CoupledSphericalSelectionRuleLinearPolarizedField_Wrapper(PyObject* py_self_, const CoupledSphericalSelectionRuleLinearPolarizedField& p0):
        CoupledSphericalSelectionRuleLinearPolarizedField(p0), py_self(py_self_) {}

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
    CoupledSphericalSelectionRuleDiagonal_Wrapper(PyObject* py_self_, const CoupledSphericalSelectionRuleDiagonal& p0):
        CoupledSphericalSelectionRuleDiagonal(p0), py_self(py_self_) {}

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

struct CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper: CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular
{
    CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper(PyObject* py_self_, const CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular& p0):
        CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular(p0), py_self(py_self_) {}

    CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper(PyObject* py_self_):
        CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular(), py_self(py_self_) {}

    bool SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular::SelectionRule(p0, p1);
    }

    PyObject* py_self;
};

struct CoupledSphericalSelectionRuleLinearPolarizedFieldAngle_Wrapper: CoupledSphericalSelectionRuleLinearPolarizedFieldAngle
{
    CoupledSphericalSelectionRuleLinearPolarizedFieldAngle_Wrapper(PyObject* py_self_, const CoupledSphericalSelectionRuleLinearPolarizedFieldAngle& p0):
        CoupledSphericalSelectionRuleLinearPolarizedFieldAngle(p0), py_self(py_self_) {}

    CoupledSphericalSelectionRuleLinearPolarizedFieldAngle_Wrapper(PyObject* py_self_):
        CoupledSphericalSelectionRuleLinearPolarizedFieldAngle(), py_self(py_self_) {}

    bool SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const CoupledSpherical::CoupledIndex& p0, const CoupledSpherical::CoupledIndex& p1) {
        return CoupledSphericalSelectionRuleLinearPolarizedFieldAngle::SelectionRule(p0, p1);
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

    class_< CoupledSphericalSelectionRuleR12, bases< CoupledSphericalSelectionRule > , CoupledSphericalSelectionRuleR12_Wrapper >("CoupledSphericalSelectionRuleR12", init<  >())
        .def(init< const CoupledSphericalSelectionRuleR12& >())
        .def(init< int >())
        .def_readwrite("cg", &CoupledSphericalSelectionRuleR12::cg)
        .def("SelectionRule", (bool (CoupledSphericalSelectionRuleR12::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&) )&CoupledSphericalSelectionRuleR12::SelectionRule, (bool (CoupledSphericalSelectionRuleR12_Wrapper::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&))&CoupledSphericalSelectionRuleR12_Wrapper::default_SelectionRule)
    ;

    class_< CoupledSphericalSelectionRuleR12Old, bases< CoupledSphericalSelectionRule > , CoupledSphericalSelectionRuleR12Old_Wrapper >("CoupledSphericalSelectionRuleR12Old", init<  >())
        .def(init< const CoupledSphericalSelectionRuleR12Old& >())
        .def_readwrite("cg", &CoupledSphericalSelectionRuleR12Old::cg)
        .def("SelectionRule", (bool (CoupledSphericalSelectionRuleR12Old::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&) )&CoupledSphericalSelectionRuleR12Old::SelectionRule, (bool (CoupledSphericalSelectionRuleR12Old_Wrapper::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&))&CoupledSphericalSelectionRuleR12Old_Wrapper::default_SelectionRule)
    ;

    class_< CoupledSphericalSelectionRuleLinearPolarizedField, bases< CoupledSphericalSelectionRule > , CoupledSphericalSelectionRuleLinearPolarizedField_Wrapper >("CoupledSphericalSelectionRuleLinearPolarizedField", init<  >())
        .def(init< const CoupledSphericalSelectionRuleLinearPolarizedField& >())
        .def("SelectionRule", (bool (CoupledSphericalSelectionRuleLinearPolarizedField::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&) )&CoupledSphericalSelectionRuleLinearPolarizedField::SelectionRule, (bool (CoupledSphericalSelectionRuleLinearPolarizedField_Wrapper::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&))&CoupledSphericalSelectionRuleLinearPolarizedField_Wrapper::default_SelectionRule)
    ;

    class_< CoupledSphericalSelectionRuleDiagonal, bases< CoupledSphericalSelectionRule > , CoupledSphericalSelectionRuleDiagonal_Wrapper >("CoupledSphericalSelectionRuleDiagonal", init<  >())
        .def(init< const CoupledSphericalSelectionRuleDiagonal& >())
        .def("SelectionRule", (bool (CoupledSphericalSelectionRuleDiagonal::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&) )&CoupledSphericalSelectionRuleDiagonal::SelectionRule, (bool (CoupledSphericalSelectionRuleDiagonal_Wrapper::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&))&CoupledSphericalSelectionRuleDiagonal_Wrapper::default_SelectionRule)
    ;

    class_< CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular, bases< CoupledSphericalSelectionRule > , CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper >("CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular", init<  >())
        .def(init< const CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular& >())
        .def("SelectionRule", (bool (CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&) )&CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular::SelectionRule, (bool (CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&))&CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper::default_SelectionRule)
    ;

    class_< CoupledSphericalSelectionRuleLinearPolarizedFieldAngle, bases< CoupledSphericalSelectionRule > , CoupledSphericalSelectionRuleLinearPolarizedFieldAngle_Wrapper >("CoupledSphericalSelectionRuleLinearPolarizedFieldAngle", init<  >())
        .def(init< const CoupledSphericalSelectionRuleLinearPolarizedFieldAngle& >())
        .def("SelectionRule", (bool (CoupledSphericalSelectionRuleLinearPolarizedFieldAngle::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&) )&CoupledSphericalSelectionRuleLinearPolarizedFieldAngle::SelectionRule, (bool (CoupledSphericalSelectionRuleLinearPolarizedFieldAngle_Wrapper::*)(const CoupledSpherical::CoupledIndex&, const CoupledSpherical::CoupledIndex&))&CoupledSphericalSelectionRuleLinearPolarizedFieldAngle_Wrapper::default_SelectionRule)
    ;

}


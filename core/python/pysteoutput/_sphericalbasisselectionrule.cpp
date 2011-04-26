
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <tensorpotential/sphericalbasisselectionrule.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

struct SphericalBasisSelectionRule_Wrapper: SphericalBasisSelectionRule
{
    SphericalBasisSelectionRule_Wrapper(PyObject* py_self_, const SphericalBasisSelectionRule& p0):
        SphericalBasisSelectionRule(p0), py_self(py_self_) {}

    SphericalBasisSelectionRule_Wrapper(PyObject* py_self_):
        SphericalBasisSelectionRule(), py_self(py_self_) {}

    bool SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    PyObject* py_self;
};

struct SphericalBasisSelectionRuleLinearPolarizedField_Wrapper: SphericalBasisSelectionRuleLinearPolarizedField
{
    SphericalBasisSelectionRuleLinearPolarizedField_Wrapper(PyObject* py_self_, const SphericalBasisSelectionRuleLinearPolarizedField& p0):
        SphericalBasisSelectionRuleLinearPolarizedField(p0), py_self(py_self_) {}

    SphericalBasisSelectionRuleLinearPolarizedField_Wrapper(PyObject* py_self_):
        SphericalBasisSelectionRuleLinearPolarizedField(), py_self(py_self_) {}

    bool SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return SphericalBasisSelectionRuleLinearPolarizedField::SelectionRule(p0, p1);
    }

    PyObject* py_self;
};

struct SphericalBasisSelectionRuleDiagonal_Wrapper: SphericalBasisSelectionRuleDiagonal
{
    SphericalBasisSelectionRuleDiagonal_Wrapper(PyObject* py_self_, const SphericalBasisSelectionRuleDiagonal& p0):
        SphericalBasisSelectionRuleDiagonal(p0), py_self(py_self_) {}

    SphericalBasisSelectionRuleDiagonal_Wrapper(PyObject* py_self_):
        SphericalBasisSelectionRuleDiagonal(), py_self(py_self_) {}

    bool SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return SphericalBasisSelectionRuleDiagonal::SelectionRule(p0, p1);
    }

    PyObject* py_self;
};

struct SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper: SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular
{
    SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper(PyObject* py_self_, const SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular& p0):
        SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular(p0), py_self(py_self_) {}

    SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper(PyObject* py_self_):
        SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular(), py_self(py_self_) {}

    bool SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular::SelectionRule(p0, p1);
    }

    PyObject* py_self;
};

struct SphericalBasisSelectionRuleLinearPolarizedFieldAngle_Wrapper: SphericalBasisSelectionRuleLinearPolarizedFieldAngle
{
    SphericalBasisSelectionRuleLinearPolarizedFieldAngle_Wrapper(PyObject* py_self_, const SphericalBasisSelectionRuleLinearPolarizedFieldAngle& p0):
        SphericalBasisSelectionRuleLinearPolarizedFieldAngle(p0), py_self(py_self_) {}

    SphericalBasisSelectionRuleLinearPolarizedFieldAngle_Wrapper(PyObject* py_self_):
        SphericalBasisSelectionRuleLinearPolarizedFieldAngle(), py_self(py_self_) {}

    bool SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return SphericalBasisSelectionRuleLinearPolarizedFieldAngle::SelectionRule(p0, p1);
    }

    PyObject* py_self;
};

struct SphericalBasisSelectionRuleDiatomicCoulomb_Wrapper: SphericalBasisSelectionRuleDiatomicCoulomb
{
    SphericalBasisSelectionRuleDiatomicCoulomb_Wrapper(PyObject* py_self_, const SphericalBasisSelectionRuleDiatomicCoulomb& p0):
        SphericalBasisSelectionRuleDiatomicCoulomb(p0), py_self(py_self_) {}

    SphericalBasisSelectionRuleDiatomicCoulomb_Wrapper(PyObject* py_self_):
        SphericalBasisSelectionRuleDiatomicCoulomb(), py_self(py_self_) {}

    SphericalBasisSelectionRuleDiatomicCoulomb_Wrapper(PyObject* py_self_, int p0):
        SphericalBasisSelectionRuleDiatomicCoulomb(p0), py_self(py_self_) {}

    bool SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return call_method< bool >(py_self, "SelectionRule", p0, p1);
    }

    bool default_SelectionRule(const SphericalBasis::LmIndex& p0, const SphericalBasis::LmIndex& p1) {
        return SphericalBasisSelectionRuleDiatomicCoulomb::SelectionRule(p0, p1);
    }

    PyObject* py_self;
};


}// namespace 


// Module ======================================================================
void Export_python_sphericalbasisselectionrule()
{
    class_< SphericalBasisSelectionRule, boost::noncopyable, SphericalBasisSelectionRule_Wrapper >("SphericalBasisSelectionRule", init<  >())
        .def("SelectionRule", pure_virtual(&SphericalBasisSelectionRule::SelectionRule))
        .def("GetBasisPairs", &SphericalBasisSelectionRule::GetBasisPairs)
    ;

    class_< SphericalBasisSelectionRuleLinearPolarizedField, bases< SphericalBasisSelectionRule > , SphericalBasisSelectionRuleLinearPolarizedField_Wrapper >("SphericalBasisSelectionRuleLinearPolarizedField", init<  >())
        .def(init< const SphericalBasisSelectionRuleLinearPolarizedField& >())
        .def("SelectionRule", (bool (SphericalBasisSelectionRuleLinearPolarizedField::*)(const SphericalBasis::LmIndex&, const SphericalBasis::LmIndex&) )&SphericalBasisSelectionRuleLinearPolarizedField::SelectionRule, (bool (SphericalBasisSelectionRuleLinearPolarizedField_Wrapper::*)(const SphericalBasis::LmIndex&, const SphericalBasis::LmIndex&))&SphericalBasisSelectionRuleLinearPolarizedField_Wrapper::default_SelectionRule)
    ;

    class_< SphericalBasisSelectionRuleDiagonal, bases< SphericalBasisSelectionRule > , SphericalBasisSelectionRuleDiagonal_Wrapper >("SphericalBasisSelectionRuleDiagonal", init<  >())
        .def(init< const SphericalBasisSelectionRuleDiagonal& >())
        .def("SelectionRule", (bool (SphericalBasisSelectionRuleDiagonal::*)(const SphericalBasis::LmIndex&, const SphericalBasis::LmIndex&) )&SphericalBasisSelectionRuleDiagonal::SelectionRule, (bool (SphericalBasisSelectionRuleDiagonal_Wrapper::*)(const SphericalBasis::LmIndex&, const SphericalBasis::LmIndex&))&SphericalBasisSelectionRuleDiagonal_Wrapper::default_SelectionRule)
    ;

    class_< SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular, bases< SphericalBasisSelectionRule > , SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper >("SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular", init<  >())
        .def(init< const SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular& >())
        .def("SelectionRule", (bool (SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular::*)(const SphericalBasis::LmIndex&, const SphericalBasis::LmIndex&) )&SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular::SelectionRule, (bool (SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper::*)(const SphericalBasis::LmIndex&, const SphericalBasis::LmIndex&))&SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular_Wrapper::default_SelectionRule)
    ;

    class_< SphericalBasisSelectionRuleLinearPolarizedFieldAngle, bases< SphericalBasisSelectionRule > , SphericalBasisSelectionRuleLinearPolarizedFieldAngle_Wrapper >("SphericalBasisSelectionRuleLinearPolarizedFieldAngle", init<  >())
        .def(init< const SphericalBasisSelectionRuleLinearPolarizedFieldAngle& >())
        .def("SelectionRule", (bool (SphericalBasisSelectionRuleLinearPolarizedFieldAngle::*)(const SphericalBasis::LmIndex&, const SphericalBasis::LmIndex&) )&SphericalBasisSelectionRuleLinearPolarizedFieldAngle::SelectionRule, (bool (SphericalBasisSelectionRuleLinearPolarizedFieldAngle_Wrapper::*)(const SphericalBasis::LmIndex&, const SphericalBasis::LmIndex&))&SphericalBasisSelectionRuleLinearPolarizedFieldAngle_Wrapper::default_SelectionRule)
    ;

    class_< SphericalBasisSelectionRuleDiatomicCoulomb, bases< SphericalBasisSelectionRule > , SphericalBasisSelectionRuleDiatomicCoulomb_Wrapper >("SphericalBasisSelectionRuleDiatomicCoulomb", init<  >())
        .def(init< const SphericalBasisSelectionRuleDiatomicCoulomb& >())
        .def(init< int >())
        .def_readwrite("cg", &SphericalBasisSelectionRuleDiatomicCoulomb::cg)
        .def("SelectionRule", (bool (SphericalBasisSelectionRuleDiatomicCoulomb::*)(const SphericalBasis::LmIndex&, const SphericalBasis::LmIndex&) )&SphericalBasisSelectionRuleDiatomicCoulomb::SelectionRule, (bool (SphericalBasisSelectionRuleDiatomicCoulomb_Wrapper::*)(const SphericalBasis::LmIndex&, const SphericalBasis::LmIndex&))&SphericalBasisSelectionRuleDiatomicCoulomb_Wrapper::default_SelectionRule)
        .def("Coefficient", &SphericalBasisSelectionRuleDiatomicCoulomb::Coefficient)
        .def("MultipoleCoeff", &SphericalBasisSelectionRuleDiatomicCoulomb::MultipoleCoeff)
        .def("CondonShortleyPhase", &SphericalBasisSelectionRuleDiatomicCoulomb::CondonShortleyPhase)
        .def("kronecker", &SphericalBasisSelectionRuleDiatomicCoulomb::kronecker)
        .staticmethod("Coefficient")
        .staticmethod("MultipoleCoeff")
        .staticmethod("CondonShortleyPhase")
        .staticmethod("kronecker")
    ;

}


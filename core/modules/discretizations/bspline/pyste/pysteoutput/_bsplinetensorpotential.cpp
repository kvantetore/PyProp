
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <tensorpotential/bsplinesolver.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
#include "tensorpotential/basis_bspline.h"


// Module ======================================================================
void Export_core_modules_discretizations_bspline_pyste_bsplinetensorpotential()
{
    class_< BSpline::BSplineSolver<1>, boost::noncopyable >("BSpline_BSplineSolver_1", init<  >())
        .def_readwrite("MatrixData", &BSpline::BSplineSolver<1>::MatrixData)
        .def_readwrite("PotentialData", &BSpline::BSplineSolver<1>::PotentialData)
        .def_readwrite("PivotData", &BSpline::BSplineSolver<1>::PivotData)
        .def("AddTensorPotential", &BSpline::BSplineSolver<1>::AddTensorPotential)
        .def("Setup", &BSpline::BSplineSolver<1>::Setup)
        .def("Solve", &BSpline::BSplineSolver<1>::Solve)
    ;

    class_< BSpline::BSplineSolver<2>, boost::noncopyable >("BSpline_BSplineSolver_2", init<  >())
        .def_readwrite("MatrixData", &BSpline::BSplineSolver<2>::MatrixData)
        .def_readwrite("PotentialData", &BSpline::BSplineSolver<2>::PotentialData)
        .def_readwrite("PivotData", &BSpline::BSplineSolver<2>::PivotData)
        .def("AddTensorPotential", &BSpline::BSplineSolver<2>::AddTensorPotential)
        .def("Setup", &BSpline::BSplineSolver<2>::Setup)
        .def("Solve", &BSpline::BSplineSolver<2>::Solve)
    ;

    class_< BSpline::BSplineSolver<3>, boost::noncopyable >("BSpline_BSplineSolver_3", init<  >())
        .def_readwrite("MatrixData", &BSpline::BSplineSolver<3>::MatrixData)
        .def_readwrite("PotentialData", &BSpline::BSplineSolver<3>::PotentialData)
        .def_readwrite("PivotData", &BSpline::BSplineSolver<3>::PivotData)
        .def("AddTensorPotential", &BSpline::BSplineSolver<3>::AddTensorPotential)
        .def("Setup", &BSpline::BSplineSolver<3>::Setup)
        .def("Solve", &BSpline::BSplineSolver<3>::Solve)
    ;

    class_< BSpline::BSplineSolver<4>, boost::noncopyable >("BSpline_BSplineSolver_4", init<  >())
        .def_readwrite("MatrixData", &BSpline::BSplineSolver<4>::MatrixData)
        .def_readwrite("PotentialData", &BSpline::BSplineSolver<4>::PotentialData)
        .def_readwrite("PivotData", &BSpline::BSplineSolver<4>::PivotData)
        .def("AddTensorPotential", &BSpline::BSplineSolver<4>::AddTensorPotential)
        .def("Setup", &BSpline::BSplineSolver<4>::Setup)
        .def("Solve", &BSpline::BSplineSolver<4>::Solve)
    ;

def("RepresentPotentialInBasisBSpline", RepresentPotentialInBasisBSpline<cplx, 1>);
def("RepresentPotentialInBasisBSpline", RepresentPotentialInBasisBSpline<cplx, 2>);
def("RepresentPotentialInBasisBSpline", RepresentPotentialInBasisBSpline<cplx, 3>);
}


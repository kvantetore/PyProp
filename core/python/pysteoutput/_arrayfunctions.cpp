
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <python/arrayfunctions.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_arrayfunctions()
{
def("SetWavefunctionFromGridFunction_1", SetWavefunctionFromGridFunction<1>);
def("SetWavefunctionFromGridFunction_2", SetWavefunctionFromGridFunction<2>);
def("SetWavefunctionFromGridFunction_3", SetWavefunctionFromGridFunction<3>);
def("SetWavefunctionFromGridFunction_4", SetWavefunctionFromGridFunction<4>);
def("SetPotentialFromGridFunction_1", SetPotentialFromGridFunction<1>);
def("SetPotentialFromGridFunction_2", SetPotentialFromGridFunction<2>);
def("SetPotentialFromGridFunction_3", SetPotentialFromGridFunction<3>);
def("SetPotentialFromGridFunction_4", SetPotentialFromGridFunction<4>);
def("OuterProduct_1", PythonOuterProduct<cplx, 1>);
def("OuterProduct_2", PythonOuterProduct<cplx, 2>);
def("OuterProduct_3", PythonOuterProduct<cplx, 3>);
def("OuterProduct_4", PythonOuterProduct<cplx, 4>);
def("OuterProduct_1", PythonOuterProduct<double, 1>);
def("OuterProduct_2", PythonOuterProduct<double, 2>);
def("OuterProduct_3", PythonOuterProduct<double, 3>);
def("OuterProduct_4", PythonOuterProduct<double, 4>);
}


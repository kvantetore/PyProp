
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <python/orthopolfunctions.cpp>
#include <utility/gamma.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_orthopolfunctions()
{
def("LaguerreQuad", PythonLaguerreQuad);
def("LaguerreMatrix", PythonLaguerreMatrix);
def("HermiteQuad", PythonHermiteQuad);
def("HermiteMatrix", PythonHermiteMatrix);
def("ScaledGridAndWeights", PythonScaledGridAndWeights);
def("OrthoPolSetup", PythonOrthoPolSetup);
def("Gamma", Gamma);
}


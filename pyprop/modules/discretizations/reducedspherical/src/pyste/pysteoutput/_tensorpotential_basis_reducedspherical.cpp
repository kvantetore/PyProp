
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
#include "src/tensorpotential/basis_reducedspherical.h"


// Module ======================================================================
void Export_pyprop_modules_discretizations_reducedspherical_src_pyste_tensorpotential_basis_reducedspherical()
{
def("RepresentPotentialInBasisReducedSphericalHarmonic", RepresentPotentialInBasisReducedSphericalHarmonic<cplx, 1>);
def("RepresentPotentialInBasisReducedSphericalHarmonic", RepresentPotentialInBasisReducedSphericalHarmonic<cplx, 2>);
def("RepresentPotentialInBasisReducedSphericalHarmonic", RepresentPotentialInBasisReducedSphericalHarmonic<cplx, 3>);
}


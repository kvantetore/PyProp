
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
#include "tensorpotential/basis_reducedspherical.h"
#include "tensorpotential/basis_finitedifference.h"


// Module ======================================================================
void Export_python_tensorpotential_basis()
{
def("RepresentPotentialInBasisReducedSphericalHarmonic", RepresentPotentialInBasisReducedSphericalHarmonic<cplx, 1>);
def("RepresentPotentialInBasisReducedSphericalHarmonic", RepresentPotentialInBasisReducedSphericalHarmonic<cplx, 2>);
def("RepresentPotentialInBasisReducedSphericalHarmonic", RepresentPotentialInBasisReducedSphericalHarmonic<cplx, 3>);
def("RepresentPotentialInBasisFiniteDifference", RepresentPotentialInBasisFiniteDifference<cplx, 1>);
def("RepresentPotentialInBasisFiniteDifference", RepresentPotentialInBasisFiniteDifference<cplx, 2>);
def("RepresentPotentialInBasisFiniteDifference", RepresentPotentialInBasisFiniteDifference<cplx, 3>);
}


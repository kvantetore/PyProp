
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
#include "tensorpotential/basis_finitedifference.h"


// Module ======================================================================
void Export_core_core_pyste_tensorpotential_basis()
{
def("RepresentPotentialInBasisFiniteDifference", RepresentPotentialInBasisFiniteDifference<cplx, 1>);
def("RepresentPotentialInBasisFiniteDifference", RepresentPotentialInBasisFiniteDifference<cplx, 2>);
def("RepresentPotentialInBasisFiniteDifference", RepresentPotentialInBasisFiniteDifference<cplx, 3>);
}


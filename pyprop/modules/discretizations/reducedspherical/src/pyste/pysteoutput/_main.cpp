// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_discretizations_reducedspherical_src_pyste_reducedsphericalkineticpotential();
void Export_pyprop_modules_discretizations_reducedspherical_src_pyste_reducedsphericalrepresentation();
void Export_pyprop_modules_discretizations_reducedspherical_src_pyste_reducedsphericaltransform();
void Export_pyprop_modules_discretizations_reducedspherical_src_pyste_tensorpotential_basis_reducedspherical();

// Module ======================================================================
BOOST_PYTHON_MODULE(libreducedspherical)
{
    Export_pyprop_modules_discretizations_reducedspherical_src_pyste_reducedsphericalkineticpotential();
    Export_pyprop_modules_discretizations_reducedspherical_src_pyste_reducedsphericalrepresentation();
    Export_pyprop_modules_discretizations_reducedspherical_src_pyste_reducedsphericaltransform();
    Export_pyprop_modules_discretizations_reducedspherical_src_pyste_tensorpotential_basis_reducedspherical();
}

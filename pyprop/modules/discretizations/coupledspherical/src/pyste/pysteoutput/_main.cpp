// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_discretizations_coupledspherical_src_pyste_coupledsphericalharmonicrepresentation();
void Export_pyprop_modules_discretizations_coupledspherical_src_pyste_coupledsphericalselectionrule();

// Module ======================================================================
BOOST_PYTHON_MODULE(libcoupledspherical)
{
    Export_pyprop_modules_discretizations_coupledspherical_src_pyste_coupledsphericalharmonicrepresentation();
    Export_pyprop_modules_discretizations_coupledspherical_src_pyste_coupledsphericalselectionrule();
}

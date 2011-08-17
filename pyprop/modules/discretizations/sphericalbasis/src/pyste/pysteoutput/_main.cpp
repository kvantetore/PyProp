// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_discretizations_sphericalbasis_src_pyste_sphericalbasisselectionrule();
void Export_pyprop_modules_discretizations_sphericalbasis_src_pyste_sphericalharmonicbasisrepresentation();

// Module ======================================================================
BOOST_PYTHON_MODULE(libsphericalbasis)
{
    Export_pyprop_modules_discretizations_sphericalbasis_src_pyste_sphericalbasisselectionrule();
    Export_pyprop_modules_discretizations_sphericalbasis_src_pyste_sphericalharmonicbasisrepresentation();
}

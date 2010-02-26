// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_discretizations_bspline_src_pyste_bsplinemain();
void Export_pyprop_modules_discretizations_bspline_src_pyste_bsplinetensorpotential();
void Export_pyprop_modules_discretizations_bspline_src_pyste_bsplinerepresentation();
void Export_pyprop_modules_discretizations_bspline_src_pyste_bsplinetransform();

// Module ======================================================================
BOOST_PYTHON_MODULE(libbspline)
{
    Export_pyprop_modules_discretizations_bspline_src_pyste_bsplinemain();
    Export_pyprop_modules_discretizations_bspline_src_pyste_bsplinetensorpotential();
    Export_pyprop_modules_discretizations_bspline_src_pyste_bsplinerepresentation();
    Export_pyprop_modules_discretizations_bspline_src_pyste_bsplinetransform();
}

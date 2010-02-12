// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_core_modules_discretizations_bspline_pyste_bsplinemain();
void Export_core_modules_discretizations_bspline_pyste_bsplinetensorpotential();
void Export_core_modules_discretizations_bspline_pyste_bsplinerepresentation();
void Export_core_modules_discretizations_bspline_pyste_bsplinetransform();

// Module ======================================================================
BOOST_PYTHON_MODULE(libbspline)
{
    Export_core_modules_discretizations_bspline_pyste_bsplinemain();
    Export_core_modules_discretizations_bspline_pyste_bsplinetensorpotential();
    Export_core_modules_discretizations_bspline_pyste_bsplinerepresentation();
    Export_core_modules_discretizations_bspline_pyste_bsplinetransform();
}

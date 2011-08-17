// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_discretizations_finitedifference_src_pyste_finitedifferencehelper();

// Module ======================================================================
BOOST_PYTHON_MODULE(libfinitedifference)
{
    Export_pyprop_modules_discretizations_finitedifference_src_pyste_finitedifferencehelper();
}

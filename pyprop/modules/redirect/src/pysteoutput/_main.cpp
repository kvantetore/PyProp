// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_redirect_src_wrapper();

// Module ======================================================================
BOOST_PYTHON_MODULE(libredirect)
{
    Export_pyprop_modules_redirect_src_wrapper();
}

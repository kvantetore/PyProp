// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_core_modules_redirect_wrapper();

// Module ======================================================================
BOOST_PYTHON_MODULE(libredirect)
{
    Export_core_modules_redirect_wrapper();
}

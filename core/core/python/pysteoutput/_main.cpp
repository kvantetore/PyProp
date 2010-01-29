// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_core_core_python_wavefunction();

// Module ======================================================================
BOOST_PYTHON_MODULE(libcore)
{
    Export_core_core_python_wavefunction();
}

// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_potentials_tensorpotentialbase_src_potential();

// Module ======================================================================
BOOST_PYTHON_MODULE(libtensorpotentialbase)
{
    Export_pyprop_modules_potentials_tensorpotentialbase_src_potential();
}

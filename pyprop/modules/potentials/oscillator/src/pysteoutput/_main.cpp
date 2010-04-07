// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_potentials_oscillator_src_potential();

// Module ======================================================================
BOOST_PYTHON_MODULE(liboscillatorpotential)
{
    Export_pyprop_modules_potentials_oscillator_src_potential();
}

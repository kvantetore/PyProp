// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_core_core_python_distributedmodel();
void Export_core_core_python_distributedoverlapmatrix();
void Export_core_core_python_wavefunction();
void Export_core_core_python_arrayfunctions();
void Export_core_core_python_blitzblas();
void Export_core_core_python_overlapmatrix();
void Export_core_core_python_blitzwrapper();
void Export_core_core_python_papi();
void Export_core_core_python_representation();
void Export_core_core_python_combinedrepresentation();
void Export_core_core_python_configurationwrapper();
void Export_core_core_python_staticpotential();
void Export_core_core_python_tensorpotential_basis();
void Export_core_core_python_customgridrepresentation();
void Export_core_core_python_tensorpotential();
void Export_core_core_python_databuffer();
void Export_core_core_python_vectorrepresentation();

// Module ======================================================================
BOOST_PYTHON_MODULE(libcore)
{
    Export_core_core_python_distributedmodel();
    Export_core_core_python_distributedoverlapmatrix();
    Export_core_core_python_wavefunction();
    Export_core_core_python_arrayfunctions();
    Export_core_core_python_blitzblas();
    Export_core_core_python_overlapmatrix();
    Export_core_core_python_blitzwrapper();
    Export_core_core_python_papi();
    Export_core_core_python_representation();
    Export_core_core_python_combinedrepresentation();
    Export_core_core_python_configurationwrapper();
    Export_core_core_python_staticpotential();
    Export_core_core_python_tensorpotential_basis();
    Export_core_core_python_customgridrepresentation();
    Export_core_core_python_tensorpotential();
    Export_core_core_python_databuffer();
    Export_core_core_python_vectorrepresentation();
}

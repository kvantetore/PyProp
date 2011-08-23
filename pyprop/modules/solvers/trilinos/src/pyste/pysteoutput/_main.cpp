// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_solvers_trilinos_src_pyste_wrapper();

// Module ======================================================================
BOOST_PYTHON_MODULE(libtrilinos)
{
    Export_pyprop_modules_solvers_trilinos_src_pyste_wrapper();
}

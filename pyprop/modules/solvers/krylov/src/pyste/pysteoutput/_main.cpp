// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_solvers_krylov_src_pyste_gmres();
void Export_pyprop_modules_solvers_krylov_src_pyste_pamp();
void Export_pyprop_modules_solvers_krylov_src_pyste_piram();

// Module ======================================================================
BOOST_PYTHON_MODULE(libkrylov)
{
    Export_pyprop_modules_solvers_krylov_src_pyste_gmres();
    Export_pyprop_modules_solvers_krylov_src_pyste_pamp();
    Export_pyprop_modules_solvers_krylov_src_pyste_piram();
}

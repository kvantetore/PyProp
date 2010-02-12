// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_core_modules_discretizations_fourier_pyste_cartesianabsorber();
void Export_core_modules_discretizations_fourier_pyste_cartesianfouriertransform();
void Export_core_modules_discretizations_fourier_pyste_cartesiankineticpotential();
void Export_core_modules_discretizations_fourier_pyste_cartesianrange();

// Module ======================================================================
BOOST_PYTHON_MODULE(libfourier)
{
    Export_core_modules_discretizations_fourier_pyste_cartesianabsorber();
    Export_core_modules_discretizations_fourier_pyste_cartesianfouriertransform();
    Export_core_modules_discretizations_fourier_pyste_cartesiankineticpotential();
    Export_core_modules_discretizations_fourier_pyste_cartesianrange();
}

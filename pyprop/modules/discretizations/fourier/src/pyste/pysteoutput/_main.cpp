// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_pyprop_modules_discretizations_fourier_src_pyste_radialtransform();
void Export_pyprop_modules_discretizations_fourier_src_pyste_cartesianabsorber();
void Export_pyprop_modules_discretizations_fourier_src_pyste_cartesianfouriertransform();
void Export_pyprop_modules_discretizations_fourier_src_pyste_cartesiankineticpotential();
void Export_pyprop_modules_discretizations_fourier_src_pyste_cartesianrange();

// Module ======================================================================
BOOST_PYTHON_MODULE(libfourier)
{
    Export_pyprop_modules_discretizations_fourier_src_pyste_radialtransform();
    Export_pyprop_modules_discretizations_fourier_src_pyste_cartesianabsorber();
    Export_pyprop_modules_discretizations_fourier_src_pyste_cartesianfouriertransform();
    Export_pyprop_modules_discretizations_fourier_src_pyste_cartesiankineticpotential();
    Export_pyprop_modules_discretizations_fourier_src_pyste_cartesianrange();
}

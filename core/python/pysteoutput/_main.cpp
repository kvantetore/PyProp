// Include =====================================================================
#include <boost/python/module.hpp>

// Exports =====================================================================
void Export_python_cranknicholson();
void Export_python_representation();
void Export_python_combinedrepresentation();
void Export_python_orthopolrepresentation();
void Export_python_orthopoltransform();
void Export_python_orthopolpotentials();
void Export_python_reducedsphericalrepresentation();
void Export_python_reducedsphericaltransform();
void Export_python_vectorrepresentation();
void Export_python_arrayfunctions();
void Export_python_orthopolfunctions();
void Export_python_blitzwrapper();
void Export_python_cartesianabsorber();
void Export_python_cartesianfouriertransform();
void Export_python_cartesiankineticpotential();
void Export_python_cartesianrange();
void Export_python_configurationwrapper();
void Export_python_distributedmodel();
void Export_python_examplepotential();
void Export_python_exponentialfinitediff();
void Export_python_sphericalabsorber();
void Export_python_sphericalkineticpotential();
void Export_python_sphericalrange();
void Export_python_sphericaltransform();
void Export_python_staticpotential();
void Export_python_polarpotential();
void Export_python_transformedgrid();
void Export_python_reducedsphericalkineticpotential();
void Export_python_wavefunction();
void Export_python_blitzblas();
void Export_python_combinedabsorber();
void Export_python_tensorpotential();
void Export_python_customgridrepresentation();

// Module ======================================================================
BOOST_PYTHON_MODULE(libcore)
{
    Export_python_cranknicholson();
    Export_python_representation();
    Export_python_combinedrepresentation();
    Export_python_orthopolrepresentation();
    Export_python_orthopoltransform();
    Export_python_orthopolpotentials();
    Export_python_reducedsphericalrepresentation();
    Export_python_reducedsphericaltransform();
    Export_python_vectorrepresentation();
    Export_python_arrayfunctions();
    Export_python_orthopolfunctions();
    Export_python_blitzwrapper();
    Export_python_cartesianabsorber();
    Export_python_cartesianfouriertransform();
    Export_python_cartesiankineticpotential();
    Export_python_cartesianrange();
    Export_python_configurationwrapper();
    Export_python_distributedmodel();
    Export_python_examplepotential();
    Export_python_exponentialfinitediff();
    Export_python_sphericalabsorber();
    Export_python_sphericalkineticpotential();
    Export_python_sphericalrange();
    Export_python_sphericaltransform();
    Export_python_staticpotential();
    Export_python_polarpotential();
    Export_python_transformedgrid();
    Export_python_reducedsphericalkineticpotential();
    Export_python_wavefunction();
    Export_python_blitzblas();
    Export_python_combinedabsorber();
    Export_python_tensorpotential();
    Export_python_customgridrepresentation();
}

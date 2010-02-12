
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <python/configuration_wrapper.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_core_core_python_configurationwrapper()
{
create_configsection_converter();
}


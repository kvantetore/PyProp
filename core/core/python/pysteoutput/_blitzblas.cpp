
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <python/blitzblas_wrapper.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_blitzblas()
{
ExportBlitzBlas();
}


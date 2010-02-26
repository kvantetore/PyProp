
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace TensorPotential { extern void export_tensorpotentialmultiply(); } 


// Module ======================================================================
void Export_core_core_src_pyste_tensorpotential()
{
TensorPotential::export_tensorpotentialmultiply();
}


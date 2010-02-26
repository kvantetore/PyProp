
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <src/redirect_funcs.cpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
BOOST_PYTHON_OPAQUE_SPECIALIZED_TYPE_ID(python_streambuf)


// Module ======================================================================
void Export_pyprop_modules_redirect_src_wrapper()
{
    def("redirect_cout", &redirect_cout, return_value_policy< return_opaque_pointer >());
    def("restore_cout", &restore_cout);
}


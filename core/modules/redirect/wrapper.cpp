
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <redirect_funcs.cpp>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
BOOST_PYTHON_OPAQUE_SPECIALIZED_TYPE_ID(python_streambuf)


// Module ======================================================================
BOOST_PYTHON_MODULE(libredirect)
{
    def("redirect_cout", &redirect_cout, return_value_policy< return_opaque_pointer >());
    def("restore_cout", &restore_cout);
}


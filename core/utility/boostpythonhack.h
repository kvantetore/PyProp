#ifndef BOOSTPYTHONHACK_H
#define BOOSTPYTHONHACK_H

#ifndef __GCCXML__
#include <boost/python.hpp>
using namespace boost::python;
#else
//UGLY hack to make pyste parse the header files. 
struct object
{
//	int a;
};
#endif

#endif


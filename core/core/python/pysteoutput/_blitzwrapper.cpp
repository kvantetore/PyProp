
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <python/array_wrapper.cpp>
#include <python/range_wrapper.cpp>
#include <python/tinyvector_wrapper.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
void Export_python_blitzwrapper()
{
import_array();
create_array_converter<double, 1>();
create_array_converter<double, 2>();
create_array_converter<double, 3>();
create_array_converter<double, 4>();
create_array_converter<int, 1>();
create_array_converter<int, 2>();
create_array_converter<int, 3>();
create_array_converter<int, 4>();
NumPyToArray<int, 1, PyArrayTraits<long> >();
NumPyToArray<int, 2, PyArrayTraits<long> >();
NumPyToArray<int, 3, PyArrayTraits<long> >();
NumPyToArray<int, 4, PyArrayTraits<long> >();
create_array_converter<std::complex<double>, 1>();
create_array_converter<std::complex<double>, 2>();
create_array_converter<std::complex<double>, 3>();
create_array_converter<std::complex<double>, 4>();
create_tinyvector_converter<int, 1>();
create_tinyvector_converter<int, 2>();
create_tinyvector_converter<int, 3>();
create_tinyvector_converter<int, 4>();
create_tinyvector_converter<double, 1>();
create_tinyvector_converter<double, 2>();
create_tinyvector_converter<double, 3>();
create_tinyvector_converter<double, 4>();
create_range_converter();
}


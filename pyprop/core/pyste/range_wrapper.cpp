
#include <boost/python/slice.hpp>

/* Converter from blitz::Range to NumPy slice */

class RangeToSlice
{
public:
	static PyObject* convert(blitz::Range range)
	{
		using namespace blitz;
		using namespace boost::python;

		int first = range.first(fromStart);
		int last = range.last(toEnd);

		object start = (first==fromStart) ? object(_) : object(first);
		object stop = (last == toEnd) ? object(_) : object(last+1); 
	
		return incref(slice(start, stop, range.stride()).ptr());
	}
};

/*
  Call this function with the appropriate template arguments to register that particular array
*/
void create_range_converter()
{
    boost::python::to_python_converter< blitz::Range, RangeToSlice >();
    //NumPyToArray<T, N>();
}




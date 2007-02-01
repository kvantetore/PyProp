
#include <Python.h>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include <blitz/array.h>
#include <complex>

#include "typedescription.h"

/* Converter from blitz::TinyVector to NumPy array */
/* Makes a copy of the objects
*/
template<class T, int N> class TinyVectorToNumPy
{
public:
	static PyObject* convert(blitz::TinyVector<T, N> vector)
	{
		using namespace blitz;
		
		npy_intp shape[1];
		npy_intp strides[1];
		shape[0] = N;
		strides[0] = sizeof(T);
		
		T v = vector(0);
		PyArray_Descr* type_descr = type_to_descr(v);
		
		//Create the array
		PyObject* wrappedObject = PyArray_NewFromDescr(&PyArray_Type, type_descr, 1, shape, strides, 0, 0, 0);
		
		//Fill it
		PyArrayObject* arrayObject = reinterpret_cast<PyArrayObject*>(wrappedObject);
		T* data = reinterpret_cast<T*>(arrayObject->data);
		for (int i=0;i<N;i++) {
			data[i] = vector(i);
		}
		
		//std::cout << vector << std::endl;
		
		return wrappedObject;
	}
};

/* Converter from NumPy Array to blitz::TinyVector */
template<class T, int N> class NumPyToTinyVector
{
public:
	NumPyToTinyVector()
	{
		//Register this class as a python type converter
		boost::python::converter::registry::push_back(
			&convertible, 
			&convert, 
			boost::python::type_id< blitz::TinyVector<T, N> >()
		);
	}
	
	static void* convertible(PyObject* obj)
	{
		if (PyArray_Check(obj))
		{
			PyArrayObject* arr_obj = (PyArrayObject*)obj;
			if (arr_obj->nd == 1) 
			{
				if (PyArray_DIM(arr_obj, 0) == N) {
					return obj;
				}
			}
		}
		
		if (PyList_Check(obj))
		{
			int size = PyList_Size(obj);
			if (size == N)
			{
				return obj;
			} else
			{
				std::cout << "Invalid length on list: " 
				          << size << " when it should have been "
					  << N 
					  << std::endl;
			}
		}
		
		return 0;
	}
	
	static void convert(PyObject* obj, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		using namespace boost::python::converter;
		using namespace boost::python;
		using namespace blitz;

		//Create vector		
		void* storage = ((rvalue_from_python_storage< TinyVector<T, N> >*) data)->storage.bytes;
		new (storage) TinyVector<T,N>();
		data->convertible = storage;
		TinyVector<T,N>* vector = reinterpret_cast<TinyVector<T,N>*>(storage);
			
		if (PyArray_Check(obj))
		{
			PyArrayObject* arr_obj = (PyArrayObject*) obj;
			T* array_data = reinterpret_cast<T*>(arr_obj->data);

			//Copy data
			for (int i=0; i<N; i++)
			{
				(*vector)(i) = array_data[i * PyArray_ITEMSIZE(arr_obj) / sizeof(T)];
			}
		}
		else if (PyList_Check(obj))
		{
			boost::python::handle<> hndl(borrowed(obj));
			object py_obj(hndl);
			extract<list> extr(py_obj);
			if (extr.check()) {
				list list_obj = extr();
				for (int i=0; i<N; i++)
				{
					object list_item = list_obj[i];
					T item_value = extract<T>(list_item);
					(*vector)(i) = item_value;
				}
			}
			else
			{
				throw std::runtime_error("This should not happen");
			}
			
		}
		else
		{
			throw std::runtime_error("How did it come to this? Trying to convert invalid type");
		}
	}
};


/*
  Call this function with the appropriate template arguments to register that particular array
*/
template<class T, int N> void create_tinyvector_converter()
{
    boost::python::to_python_converter<blitz::TinyVector<T, N>, TinyVectorToNumPy<T, N> >();
    NumPyToTinyVector<T, N>();
}


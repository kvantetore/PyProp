
#include <Python.h>
#include <numpy/arrayobject.h>
#include <blitz/array.h>
#include <complex>

#include "typedescription.h"

/* Converter from blitz::Array to NumPy array */
/* NB! only pass array which lifetime exceeds the lifetime
of this function, otherwise the NumPy array's will be freed when reference 
counter of the blitz::Array reaches zero.

It would be nice to tie up the lifespan of the array memory to both the
NumPy and blitz:: arrays, but i'm not sure how...
*/
template<class T, int N> class ArrayToNumPy
{
public:
	static PyObject* convert(blitz::Array<T, N> array)
	{
		using namespace blitz;
		
		npy_intp shape[N];
		npy_intp strides[N];
		
		T v = array(0);
		PyArray_Descr* type_descr = type_to_descr(v);
		
		//std::cout << "ref: " << array.getReferenceCount() << std::endl;
		
		for (int i=0;i<N;i++) {
			shape[i] = array.extent(i);
			strides[i] = array.stride(i) * sizeof(T);
			//std::cout << "shape: " << shape[i] << " stride: " << strides[i] << std::endl;
		}
		
		//std::cout << array << std::endl;
		
		PyObject* wrappedObject = PyArray_NewFromDescr(&PyArray_Type, type_descr, N, shape, strides, array.data(), NPY_WRITEABLE, 0);
		
		//If the Array is about to expire (ie it is not stored anywhere but in the return
		//argument of a function), we return a copy which is to be maintained by python
		if (array.getReferenceCount() <= 2) 
		{
			std::cout << "Returning python-owning copy of array" << std::endl;
			
			//Create a copy of the array (perhaps there is a more elegant way?)
			PyObject* newObject = PyArray_FromAny(wrappedObject, type_descr, 0, 0, NPY_ENSURECOPY, 0);
			
			//dispose of the old arrayObject
			Py_DECREF(wrappedObject);
			
			//return the copy of the array
			wrappedObject = newObject;
		}
		
		return wrappedObject;
	}
};

/* Converter from NumPy Array to blitz::Array */
template<class T, int N> class NumPyToArray
{
public:
	NumPyToArray()
	{
		//Register this class as a python type converter
		boost::python::converter::registry::push_back(
			&convertible, 
			&convert, 
			boost::python::type_id< blitz::Array<T, N> >()
		);
	}
	
	static void* convertible(PyObject* obj)
	{
		if (PyArray_Check(obj))
		{
			PyArrayObject* arr_obj = (PyArrayObject*)obj;
			if (arr_obj->nd == N) 
			{
				return obj;
			}
		}
		
		return 0;
	}
	
	static void convert(PyObject* obj, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		using namespace boost::python::converter;
		using namespace blitz;
		
		PyArrayObject* arr_obj = (PyArrayObject*) obj;
		
		void* storage = ((rvalue_from_python_storage< Array<T, N> >*) data)->storage.bytes;
		
		TinyVector<int,N> shape(0);
		TinyVector<int,N> strides(0);
	
		for (int i=0;i<N;i++) {
			shape[i] = arr_obj->dimensions[i];
			strides[i] = arr_obj->strides[i] / sizeof(T);
		}
	
		//Find type descriptor for C++ type
		T v = 0;
		PyArray_Descr* type_descr = type_to_descr(v);
		
		if (arr_obj->descr->type_num == type_descr->type_num) {
			new (storage) Array<T,N>((T*)arr_obj->data, shape, strides, neverDeleteData);
		}
		else
		{
			std::cout << "Invalid datatype" << std::endl;
			throw std::runtime_error("Invalid datatype");
		}
		data->convertible = storage;
	}
};

/*
  Call this function with the appropriate template arguments to register that particular array
*/
template<class T, int N> void create_array_converter()
{
    boost::python::to_python_converter<blitz::Array<T, N>, ArrayToNumPy<T, N> >();
    NumPyToArray<T, N>();
}


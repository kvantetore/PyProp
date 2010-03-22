
#include <Python.h>
#include <numpy/arrayobject.h>
#include <complex>
#include <iostream>

#include "typedescription.h"

typedef std::complex<double> cplx;
const cplx I = cplx(0.0, 1.0);

/* Converter from NumPy Array to blitz::Array */
template<class T, class Traits=PyArrayTraits<T> > class NumPyScalarToScalar
{
public:
	NumPyScalarToScalar()
	{
		//Register this class as a python type converter
		boost::python::converter::registry::push_back(
			&convertible, 
			&convert, 
			boost::python::type_id< T >()
		);
	}
	
	static void* convertible(PyObject* obj)
	{
		std::cout << "Type: " << Py_TYPE(obj) << std::endl;
		std::cout << "Type: " << &PyArray_Type << std::endl;
		if (!PyArray_Check(obj))
		{
			std::cout << "not a pyarray object" << std::endl;
			return 0;
		}

		if (!PyArray_CheckScalar(obj))
		{
			std::cout << "not a pyarray scalar" << std::endl;
			return 0;
		}

		PyArrayObject* arr_obj = (PyArrayObject*)obj;
		PyArray_Descr* from_type = arr_obj->descr;
		PyArray_Descr* to_type = Traits::GetTypeDescr();
		
		if (from_type->type_num == to_type->type_num)
		{
			return obj;
		}
		else
		{
			std::cout << "From Type " << from_type->type_num << " -> " << to_type->type_num << std::endl;
		}
	
		if (from_type->elsize == to_type->elsize)
		{
			return obj;
		}
		else
		{
			std::cout << "From elsize " << from_type->elsize  << " -> " << to_type->elsize  << std::endl;
		}

		return 0;
	}
	
	static void convert(PyObject* obj, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		using namespace boost::python::converter;
		
		PyArrayObject* arr_obj = (PyArrayObject*) obj;
		
		void* storage = ((rvalue_from_python_storage< T >*) data)->storage.bytes;
		(*(T*)storage) = (*(T*)arr_obj->data);
		//PyArray_ScalarAsCtype(obj, storage);
		data->convertible = storage;
	}
		
};

/*
  Call this function with the appropriate template arguments to register that particular array
*/
void create_array_scalar_converter()
{
    NumPyScalarToScalar<int>();
    NumPyScalarToScalar<long>();
    NumPyScalarToScalar<double>();
    NumPyScalarToScalar<cplx>();
    
}




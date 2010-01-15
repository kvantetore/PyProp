
#include <iostream>
#include <boost/python.hpp>
#include "configuration.h"

/* Converter for ConfigSection */

/* Converter from Python Section-object to C++ ConfigSection */
class ConfigSectionCppToPython
{
public:
	static PyObject* convert(const ConfigSection &section)
	{
		using namespace boost::python;
		
		object obj = *(object*)section.SectionObject;
		return obj.ptr();
	}
};

class ConfigSectionPythonToCpp
{
public:
	ConfigSectionPythonToCpp()
	{
		//Register this class as a python type converter
		boost::python::converter::registry::push_back(
			&convertible, 
			&convert, 
			boost::python::type_id< ConfigSection >()
		);
	}
	
	static void* convertible(PyObject* obj)
	{
		return obj;
	}
	
	static void convert(PyObject* obj, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		using namespace boost::python::converter;
		using namespace boost::python;
		
		void* storage = ((rvalue_from_python_storage< ConfigSection >*) data)->storage.bytes;
		new (storage) ConfigSection(obj);
		data->convertible = storage;
	}
};

/*
  Call this function to register the converter
*/
void create_configsection_converter()
{
    boost::python::to_python_converter<ConfigSection, ConfigSectionCppToPython>();
    ConfigSectionPythonToCpp();
}


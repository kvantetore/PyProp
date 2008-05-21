#include "configuration.h"
#include "common.h"

//For some reason both std::string and boost::python defines _POSIX_C_SOURCE.
//This is my little hack to ignore the warning
#include <string>
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif

#include <blitz/array.h>
#include <complex>
#include <map>

typedef std::map<std::string, std::string> ConfigMap;

/* implementations */

bool HasPythonValue(const std::string &name, const ConfigSection *conf)
{
	handle<> hndl(borrowed((PyObject*)conf->SectionObject));
	object config_obj(hndl);
	return extract<bool>(config_obj.attr("Exists")(name));
}

bool HasMapValue(const std::string &name, const ConfigSection *conf)
{
	ConfigMap *map = static_cast<ConfigMap*>(conf->SectionObject);
	return map->find(name) == map->end();
}

template<class T> T GetPythonValue(const std::string &name, const ConfigSection *conf)
{
	handle<> hndl(borrowed((PyObject*)conf->SectionObject));
	object config_obj(hndl);
	object obj = config_obj.attr(name.c_str());
	T ret = extract<T>(obj);
	
	return ret;
}

template<class T> T GetMapValue(const std::string &name, const ConfigSection *conf)
{
	ConfigMap *map = static_cast<ConfigMap*>(conf->SectionObject);
	std::string str_value = (*map)[name];
	
	//TODO:Implement fromString
	return FromString<T>(str_value);	
}


ConfigSection::ConfigSection(void* sectionObject) : 
	Mode(ConfigSectionModePython),
	SectionObject(sectionObject)
{
	Py_INCREF((PyObject*)SectionObject);
}

ConfigSection::~ConfigSection()
{
	if (SectionObject != 0) 
	{
		if (Mode == ConfigSectionModePython)
		{
			Py_DECREF((PyObject*)SectionObject);
			SectionObject = 0;
		}
	}
} 

ConfigSection::ConfigSection(const ConfigSection& other)
{
	Mode = other.Mode;
	if (Mode == ConfigSectionModePython)
	{
		//We must create a instance of object for each instance of 
		//ConfigSection in order to keep the reference count of the
		//python object.
		SectionObject = other.SectionObject;
		Py_INCREF((PyObject*)SectionObject);
	}
}

const object ConfigSection::GetPythonConfigSection() const
{
	if (Mode != ConfigSectionModePython)
	{
		throw std::runtime_error("ConfigSection Not implemented!");
	}
	const handle<> hndl(borrowed((PyObject*)this->SectionObject));
	const object config_obj(hndl);
	
	return config_obj;
}

template<class T>
T ConfigSection::Get(const std::string &name) const
{
	if (Mode == ConfigSectionModePython)
	{
		return GetPythonValue<T>(name, this);
	}
	if (Mode == ConfigSectionModeMap)
	{
		return GetMapValue<T>(name, this);
	}
	throw std::runtime_error("Invalid config section mode");
}

bool ConfigSection::HasValue(const std::string &name) const
{
	if (Mode == ConfigSectionModePython)
	{
		return HasPythonValue(name, this);
	}
	if (Mode == ConfigSectionModeMap)
	{
		return HasMapValue(name, this);
	}
	throw std::runtime_error("Invalid config section mode");
}

/*template<class T>
void ConfigSection::Set(const std::string &name, const T &value)
{
	if (Mode == ConfigSectionModePython)
	{
		throw std::runtime_error("Cannot set config value when in python mode");
	}
	if (Mode == ConfigSectionModeMap)
	{
		//TODO: Implement fromstring
		std::string str_value = ToString(value);
		GetMapObject(this)[name] = str_value;
	}
}
*/

template bool ConfigSection::Get<bool>(const std::string &name) const;
template object ConfigSection::Get<object>(const std::string &name) const;
template int ConfigSection::Get<int>(const std::string &name) const;
template double ConfigSection::Get<double>(const std::string &name) const;
template std::string ConfigSection::Get<std::string>(const std::string &name) const;

template blitz::TinyVector<int, 1> ConfigSection::Get< blitz::TinyVector<int, 1> >(const std::string &name) const;
template blitz::TinyVector<int, 2> ConfigSection::Get< blitz::TinyVector<int, 2> >(const std::string &name) const;
template blitz::TinyVector<int, 3> ConfigSection::Get< blitz::TinyVector<int, 3> >(const std::string &name) const;
template blitz::TinyVector<int, 4> ConfigSection::Get< blitz::TinyVector<int, 4> >(const std::string &name) const;
template blitz::TinyVector<int, 5> ConfigSection::Get< blitz::TinyVector<int, 5> >(const std::string &name) const;
template blitz::TinyVector<int, 6> ConfigSection::Get< blitz::TinyVector<int, 6> >(const std::string &name) const;
template blitz::TinyVector<int, 7> ConfigSection::Get< blitz::TinyVector<int, 7> >(const std::string &name) const;
template blitz::TinyVector<int, 8> ConfigSection::Get< blitz::TinyVector<int, 8> >(const std::string &name) const;
template blitz::TinyVector<int, 9> ConfigSection::Get< blitz::TinyVector<int, 9> >(const std::string &name) const;

template blitz::TinyVector<double, 1> ConfigSection::Get< blitz::TinyVector<double, 1> >(const std::string &name) const;
template blitz::TinyVector<double, 2> ConfigSection::Get< blitz::TinyVector<double, 2> >(const std::string &name) const;
template blitz::TinyVector<double, 3> ConfigSection::Get< blitz::TinyVector<double, 3> >(const std::string &name) const;
template blitz::TinyVector<double, 4> ConfigSection::Get< blitz::TinyVector<double, 4> >(const std::string &name) const;
template blitz::TinyVector<double, 5> ConfigSection::Get< blitz::TinyVector<double, 5> >(const std::string &name) const;
template blitz::TinyVector<double, 6> ConfigSection::Get< blitz::TinyVector<double, 6> >(const std::string &name) const;
template blitz::TinyVector<double, 7> ConfigSection::Get< blitz::TinyVector<double, 7> >(const std::string &name) const;
template blitz::TinyVector<double, 8> ConfigSection::Get< blitz::TinyVector<double, 8> >(const std::string &name) const;
template blitz::TinyVector<double, 9> ConfigSection::Get< blitz::TinyVector<double, 9> >(const std::string &name) const;

template blitz::Array<double, 1> ConfigSection::Get< blitz::Array<double, 1> >(const std::string &name) const;
template blitz::Array<double, 2> ConfigSection::Get< blitz::Array<double, 2> >(const std::string &name) const;
template blitz::Array<double, 3> ConfigSection::Get< blitz::Array<double, 3> >(const std::string &name) const;

template blitz::Array<cplx, 1> ConfigSection::Get< blitz::Array<cplx, 1> >(const std::string &name) const;
template blitz::Array<cplx, 2> ConfigSection::Get< blitz::Array<cplx, 2> >(const std::string &name) const;
template blitz::Array<cplx, 3> ConfigSection::Get< blitz::Array<cplx, 3> >(const std::string &name) const;

template blitz::Array<int, 1> ConfigSection::Get< blitz::Array<int, 1> >(const std::string &name) const;

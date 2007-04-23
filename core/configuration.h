#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>

enum ConfigSectionMode
{
	ConfigSectionModePython = 1,
	ConfigSectionModeMap    = 2
};

class ConfigSection
{
public:
	ConfigSectionMode Mode;
	void* SectionObject;

	ConfigSection() : Mode(ConfigSectionModeMap), SectionObject(0) {}	
	ConfigSection(void* sectionObject); 
	ConfigSection(const ConfigSection & copy); 
	virtual ~ConfigSection();

	bool HasValue(const std::string &name) const;
	template<class T> T Get(const std::string &name) const;
	template<class T> T Set(const std::string &name, const T& value);
	
	template<class T> void Get(const std::string &name, T& dest) const
	{
		dest = Get<T>(name);
	}
};

#endif

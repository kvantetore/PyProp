#ifndef PYPROP_USE_PAPI

void ExportPapi()
{
}

#else

#include <stdexcept>
#include <papi.h>

void PapiSetup()
{

	/* Initialize PAPI library */
	if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
		throw std::runtime_error("Error initializing PAPI");
	}
}

PAPI_dmem_info_t GetDynamicMemoryInfo()
{
	PAPI_dmem_info_t meminfo;
	PAPI_get_dmem_info(&meminfo);
	return meminfo;
}

void ExportPapi()
{
	using namespace boost::python;
	
	class_<PAPI_dmem_info_t>("DynamicMemoryInfo")
		.def_readonly("peak", &PAPI_dmem_info_t::peak)
		.def_readonly("size", &PAPI_dmem_info_t::size)
		.def_readonly("resident", &PAPI_dmem_info_t::resident)
		.def_readonly("high_water_mark", &PAPI_dmem_info_t::high_water_mark)
		.def_readonly("shared", &PAPI_dmem_info_t::shared)
		.def_readonly("text", &PAPI_dmem_info_t::text)
		.def_readonly("library", &PAPI_dmem_info_t::library)
		.def_readonly("heap", &PAPI_dmem_info_t::heap)
		.def_readonly("locked", &PAPI_dmem_info_t::locked)
		.def_readonly("stack", &PAPI_dmem_info_t::stack)
		.def_readonly("pagesize", &PAPI_dmem_info_t::pagesize)
		.def_readonly("pte", &PAPI_dmem_info_t::pte)
	;

	def("PapiGetDynamicMemoryInfo", GetDynamicMemoryInfo);
	def("PapiSetup", PapiSetup);
}


#endif


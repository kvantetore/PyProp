#include "angularrepresentation.h"

void AngularRepresentation::SetupRepresentation(int maxl)
{
	// For now gridtype is. Later it should be possible to
	// choose between gauss-legendre and other sets
	OmegaRange::GridType gridtype = OmegaRange::Gauss;

	//need to determine how many points given maxl and angular type.
	Range.SetupRange(gridtype, maxl);
}

void AngularRepresentation::ApplyConfigSection(const ConfigSection &config)
{
		int maxl;
		config.Get("maxl", maxl);
		SetupRepresentation(maxl);
}	


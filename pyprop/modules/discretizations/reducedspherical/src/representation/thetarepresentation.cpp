#include "thetarepresentation.h"

namespace ReducedSpherical
{

void ThetaRepresentation::SetupRepresentation(int maxl)
{
	Range.SetupRange(maxl);
}

void ThetaRepresentation::ApplyConfigSection(const ConfigSection &config)
{
		int maxl;
		config.Get("maxl", maxl);
		SetupRepresentation(maxl);
}	

} //namespace


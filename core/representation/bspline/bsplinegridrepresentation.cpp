#include "bsplinegridrepresentation.h"

namespace BSpline
{

void BSplineGridRepresentation::ApplyConfigSection(const ConfigSection &config)
{
	// Do stuff here
}

void BSplineGridRepresentation::SetupRepresentation(BSpline::Ptr p)
{
	BSplineObject = p;
	Grid = BSplineObject->GetQuadratureGridGlobal();
	Weights = BSplineObject->GetQuadratureGridGlobal();
}	

} //namespace


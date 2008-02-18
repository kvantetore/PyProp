#include "bsplinerepresentation.h"

namespace BSpline
{

/*
 *  Apply config
 */
void BSplineRepresentation::ApplyConfigSection(const ConfigSection &config)
{
	/*
	 * First we call the b-spline setup function as given
	 * in the config section. This will return a b-spline
	 * object from which all representation data is drawn.
	 */
	object initFunction = config.Get<object>("init_function");
	object pythonConf = config.GetPythonConfigSection();
	BSplineObject = extract<BSpline::Ptr>(initFunction(pythonConf));

	// Get grid/weight size from BSpline object
	int gridSize = BSplineObject->NumberOfBSplines;
	Grid.resize(gridSize);
	Grid = blitz::tensor::i;
	Weights.resize(gridSize);
	Weights = 1.0;
}

} //namespace


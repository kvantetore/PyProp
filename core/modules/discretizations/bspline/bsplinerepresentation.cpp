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
	BSpline::Ptr bsplineObject = extract<BSpline::Ptr>(initFunction(pythonConf));

	SetupRepresentation(bsplineObject);
}

void BSplineRepresentation::SetupRepresentation(BSpline::Ptr bsplineObject)
{
	BSplineObject = bsplineObject;

	// Get grid/weight size from BSpline object
	int gridSize = BSplineObject->NumberOfBSplines;
	//cout << "k = " << BSplineObject->MaxSplineOrder << " N = " << BSplineObject->NumberOfBSplines << endl;
	cout << "Got BSpline basis size " << gridSize << endl;
	cout << "Got number of BSpline integration points " << BSplineObject->GetQuadratureGridGlobal().extent(0) << endl;

	//Grid is just a range integers (0,1,2...)
	Grid.resize(gridSize);
	Grid = blitz::tensor::i;

	//Integration weights are all 1.0
	Weights.resize(gridSize);
	Weights = 1.0;

}

} //namespace


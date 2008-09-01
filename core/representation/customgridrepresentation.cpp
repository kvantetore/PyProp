#include "customgridrepresentation.h"

void CustomGridRepresentation::ApplyConfigSection(const ConfigSection &config)
{
	object gridFunction = config.Get<object>("function");
	object pythonConf = config.GetPythonConfigSection();

	blitz::Array<double, 1> grid = extract< blitz::Array<double, 1> >(gridFunction(pythonConf));
	GlobalGrid.reference(grid.copy());

	//Set up trapezoidal integration rule with zero boundary conditions
	int N = GlobalGrid.size();
	GlobalWeights.resize(N);
	for (int i=1; i<N-1; i++)
	{
		GlobalWeights(i) = 0.5 * (GlobalGrid(i+1) - GlobalGrid(i-1));
	}
	//Assume that the 0-valued boundary point is equidistant with the two previous points
	GlobalWeights(0) =  (GlobalGrid(1) - GlobalGrid(0));
	GlobalWeights(N-1) =  (GlobalGrid(N-1) - GlobalGrid(N-2));
}




#include "customgridrepresentation.h"

void CustomGridRepresentation::ApplyConfigSection(const ConfigSection &config)
{
	object gridFunction = config.Get<object>("function");
	object pythonConf = config.GetPythonConfigSection();

	blitz::Array<double, 1> grid = extract< blitz::Array<double, 1> >(gridFunction(pythonConf));
	GlobalGrid.reference(grid.copy());

	std::string quadrature = config.Get<std::string>("quadrature");
	if (quadrature == "trapez")
	{
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
	else if (quadrature == "simpson")
	{
		int N = GlobalGrid.size();
		GlobalWeights.resize(N);
		cout << "N%2 = " << (N%2) << endl;
		if (N % 2 != 1)
		{
			throw std::runtime_error("simpson requires odd number of gridpoints");
		}
		double h = GlobalGrid(1) - GlobalGrid(0);
		GlobalWeights(0) = h/3.;
		GlobalWeights(N-1) = h/3.;
		GlobalWeights(N-2) = 4*h/3.;
		for (int i=0; i<(N-3)/2; i++)
		{
			GlobalWeights(2*i + 1) = 4. * h / 3.;
			GlobalWeights(2*i + 2) = 2. * h / 3.;
		}

		cout << "GlobalWeights = " << GlobalWeights << endl;
	}
}



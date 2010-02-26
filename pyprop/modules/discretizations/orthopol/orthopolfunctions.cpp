#include <boost/python.hpp>
#include "../common.h"
#include "../transform/orthopol/orthopoltools.h"

using namespace boost::python;
using namespace blitz;
using namespace OrthoPol;

void PythonLaguerreQuad(int N, double m, Array<double, 1> X, Array<double, 1> W)
{
	LaguerreQuad(N, m, X, W);
}

void PythonLaguerreMatrix(int N, double m, Array<double, 1> X, Array<double, 2> L)
{
	LaguerreMatrix(N, m, X, L);
}

void PythonHermiteQuad(int N, Array<double, 1> X, Array<double, 1> W)
{
	HermiteQuad(N, X, W);
}

void PythonHermiteMatrix(int N, Array<double, 1> X, Array<double, 2> H)
{
	HermiteMatrix(N, X, H);
}

void PythonScaledGridAndWeights(int N, Parameter &param, Array<double, 1> X, Array<double, 1> W)
{
	ScaledGridAndWeights(N, param, X, W);
}

void PythonOrthoPolSetup(int N, Parameter &param, Array<double, 1> eigenvalues, Array<double, 2> eigenvectors, Array<double, 1> x, Array<double, 1> w)
{
	OrthoPolSetup(N, param, eigenvalues, eigenvectors, x, w);
}



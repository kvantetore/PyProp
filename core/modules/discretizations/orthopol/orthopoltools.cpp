#include "orthopoltools.h"
#include "../../utility/orthogonalpolynomial.h"
#include "../../utility/gamma.h"

namespace OrthoPol
{
	
void LaguerreQuad(int N, const double &m, Array<double, 1> &X, Array<double, 1> &W)
{
	Array<double, 1> D(N);
	Array<double, 1> E(N-1);
	D  = 2.0*tensor::i + m + 1.0;
	E  = -sqrt((tensor::i + 1) * (tensor::i + 1 + m)); 
	
	Array<double, 2> quadrature = OrthogonalPolynomial::CalculateQuadrature(N, D, E);
	X = quadrature(Range::all(), 0);
	W = quadrature(Range::all(), 1) * Gamma(1 + m) / 2.0;
}

void LaguerreMatrix(int N, const double &m, const Array<double, 1> &X, Array<double, 2> &L)
{
	Array<double,2> l(N+1, N);
	l=0;
	l(1,Range::all())=exp(-(X/2))/sqrt(Gamma(1 + m)); // 1.772453850905516 = gamma(0.5)
	
	for (int i=1; i<N; i++)
	{
		l(i+1,Range::all())=(2.0*i+m-1.0-X)*l(i,Range::all()) - sqrt((i-1.0)*(i+m-1.0))*l(i-1,Range::all());
		l(i+1,Range::all())=l(i+1,Range::all())/sqrt(i*(i+m));
	}	
	L=l(Range(1,N),Range::all());
}

void HermiteQuad(int N, Array<double, 1> &X, Array<double, 1> &W)
{
	Array<double, 1> D(N);
	Array<double, 1> E(N-1);
	D = 0;
	E = sqrt((tensor::i + 1.0)/2.0);

	Array<double, 2> quadrature = OrthogonalPolynomial::CalculateQuadrature(N, D, E);
	X = quadrature(Range::all(), 0);
	W = quadrature(Range::all(), 1) * sqrt(M_PI) / 2.0;
}

void HermiteMatrix(int N, const Array<double, 1> &X, Array<double, 2> &H)
{
	Array<double,2> h(N+1, N);
	h=0;
	h(1,Range::all())=exp(-(sqr(X)/2))/sqrt(sqrt(M_PI));

	for (int i=1; i<N; i++)
	{
		h(i+1,Range::all())=(X)*sqrt(2.0/i)*h(i,Range::all()) - sqrt((i-1.0)/i)*h(i-1,Range::all());
	}	
	H=h(Range(1,N),Range::all());
}

/*
 * Calculates grid and weights for orthopol representation
 * Input parameters:
 * 	- N:        Number of grid points
 * 	- param:    Parameters to the orthogonal polynomial
 * Output parameters:
 *  - X:        Quadrature grid points (resized to N)
 *  - W:        Quadrature weights (resized to N)
 */
void GridAndWeights(int N, const Parameter &param, Array<double, 1> &X, Array<double, 1> &W)
{
	X.resize(N);
	W.resize(N);
	if (param.Type == HermitePolynomial)
	{
		HermiteQuad(N, X, W);
		W *= exp(X*X);
	}
	else if (param.Type == LaguerrePolynomial)
	{
		//double m = -0.5;
		double m = param.HypersphericalRank/2.0 - 1.0;
		LaguerreQuad(N, m, X, W);
		W *= exp(X); //*pow(X,-m);
	}
	else 
	{ 
		cout << "Illegal parameter type!" << endl; 
	}
}

/*
 * Scales the grid x according to param.Scaling
 */
void ScaleGrid(const Parameter &param, Array<double, 1> &x)
{
	if (param.Type == HermitePolynomial)
	{
		x  = x / param.Scaling; 
	}
	else if (param.Type == LaguerrePolynomial)
	{
		x  = sqrt(x) * 2 / param.Scaling; 
		//x = x / param.Scaling;
	}
	else 
	{ 
		cout << "Illegal parameter type!" << endl; 
	}
}

void ScaledGridAndWeights(int N, const Parameter &param, Array<double, 1> &x, Array<double, 1> &w)
{
	GridAndWeights(N, param, x, w);
	ScaleGrid(param, x);
}

void OrthoPolSetup(int N, const Parameter &param, Array<double, 1> &eigenvalues, Array<double, 2> &eigenvectors, Array<double, 1> &x, Array<double, 1> &w)
{
	eigenvectors.resize(N,N);
	eigenvalues.resize(N);

	// Calculate Eigenvalues
	eigenvalues = sqr(param.Scaling) * tensor::i;

	// Find Quadrature points and weights.
	GridAndWeights(N, param, x, w);

	//Calculate the sampled functions by 
	if (param.Type== HermitePolynomial)
	{
		HermiteMatrix(N, x, eigenvectors);
	}
	else if (param.Type== LaguerrePolynomial)
	{
		//double m = -0.5;
		double m = param.HypersphericalRank/2.0 - 1.0;
		LaguerreMatrix(N, m, x, eigenvectors);
		eigenvalues /= 2.0;
	}
	else 
	{ 
		cout << "Illegal parameter type!" << endl; 
	}

	//Scale the grid
	ScaleGrid(param, x);
	
}	

}; //Namespae


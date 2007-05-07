#include "orthopoltools.h"
#include "../../utility/orthogonalpolynomial.h"

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
	W = quadrature(Range::all(), 1) * 1.772453850905516 / 2.0;
}

void LaguerreMatrix(int N, const double &m, const Array<double, 1> &X, Array<double, 2> &L)
{
	Array<double,2> l(N+1, N);
	l=0;
	l(1,Range::all())=exp(-(X/2))/sqrt(1.772453850905516); // 1.772453850905516 = gamma(0.5)
	
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
	h(1,Range::all())=exp(-(sqr(X)/2))/sqrt(sqrt(3.141592653589793));

	for (int i=1; i<N; i++)
	{
		h(i+1,Range::all())=(X)*sqrt(2.0/i)*h(i,Range::all()) - sqrt((i-1.0)/i)*h(i-1,Range::all());
	}	
	H=h(Range(1,N),Range::all());
}

void SetupLaguerre(int N, const cplx &dt, const double &cutoff, Array<cplx, 2> &P)
{
	Array<double, 1> w(N);
	Array<double, 1> X(N);
	Array<double, 2> L(N, N);
	
	// Ass. Laguerre: L^(m) (x).
	double m = -0.5;

	// Find Quadrature points and weights.
	LaguerreQuad(N, m, X, w);
	
	// Compute Laguerre Functions.
	LaguerreMatrix(N, m, X, L);
	
	// Stretch variable.
	double alpha = sqrt(max(X))*2/cutoff;
	
	// Compute Propagation Matrix, P.
	Array<cplx, 1> d(N);
	w *= exp(X);
	for (int i=0; i<N; i++)
	{
		d(i) = exp(-I*dt*sqr(alpha)*(double)i/2.0);
	}
	cplx dot = 0;
	
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			for (int k=0; k<N; k++)
			{
				dot += L(k,i)*d(k)*L(k,j);
			}
			dot *= w(j);
			P(i,j) = dot;
			dot = 0;
		}
	}
}	


void SetupHermite(int N, const cplx &dt, const double &cutoff, Array<cplx, 2> &P)
{
	
	Array<double, 1> w(N);
	Array<double, 1> X(N);
	Array<double, 2> H(N, N);
	
	// Find Quadrature points and weights.
	HermiteQuad(N, X, w);

	// Compute Hermite Functions.
	HermiteMatrix(N, X, H);
	
	// Stretch variable.
	double alpha = max(X)/cutoff;
	
	// Compute Propagation Matrix, P.
	Array<cplx, 1> d(N);
	w *= exp(sqr(X));
	for (int i=0; i<N; i++)
	{
		d(i) = exp(-I*dt*sqr(alpha)*(double)i);
	}
	cplx dot = 0;
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			dot = 0;
			for (int k=0; k<N; k++)
			{
				dot += H(k,i)*d(k)*H(k,j);
			}
			dot *= w(j);
			P(i,j) = dot;
		}
	}
}


/*
 * Calculates grid and weights for orthopol representation
 * Input parameters:
 * 	- N:        Number of grid points
 * 	- Param:    Parameters to the orthogonal polynomial
 * Output parameters:
 *  - X:        Quadrature grid points (resized to N)
 *  - W:        Quadrature weights (resized to N)
 *  - Alpha:    Grid scaling
 */
void GridAndWeights(int N, const Parameter &Param, Array<double, 1> &X, Array<double, 1> &W, double &Alpha)
{
	X.resize(N);
	W.resize(N);
	double L = Param.Cutoff;
	if (Param.Type == HermiteTransform)
	{
		HermiteQuad(N, X, W);
		Alpha = max(X) / L;
		W *= exp(X*X);
		X  = X / Alpha; // new X defined on [-L,L]
	}
	else if (Param.Type == LaguerreTransform)
	{
		double m = -0.5;
		LaguerreQuad(N, m, X, W);
		W *= exp(X)*pow(X,-m);
		Alpha = 2 * sqrt(max(X)) / L;
		X  = 2 * sqrt(X) / Alpha; // new X defined on (0,L]
	}
	else 
	{ 
		cout << "Illegal parameter type!" << endl; 
	}
}

void OrthoPolSetup(int N, const Parameter &Param, const cplx &dt, Array<cplx, 2> &P)
{
	P.resize(N,N);
	if (Param.Type== HermiteTransform)
	{
		SetupHermite(N, dt, Param.Cutoff, P);
	}
	else if (Param.Type== LaguerreTransform)
	{
		SetupLaguerre(N, dt, Param.Cutoff, P);
	}
	else 
	{ 
		cout << "Illegal parameter type!" << endl; 
	}
}	

}; //Namespae


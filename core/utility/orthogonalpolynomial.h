#include "../common.h"

namespace OrthogonalPolynomial
{

enum PolynomialType
{
	Legendre1,        //Legendre polys on (-1, 1)
	Legendre2,        //Legendre polys on (0, 1)
	Chebyshev1, 	  //Chebyshev polys of first kind
	Chebyshev2,       //Chebyshev polys of second kind
	PolynomialTypeEnd
};

/*
 * Calculates the weights and grid points for the N-point gaussian quadrature.
 * diagonal and subDiagonal are the specification of the recursion coefficients
 * to the corresponding orthogonal polynomial
 * It returns a [Nx2] array quad, where
 * quad(:, 0) == grid points
 * quaD(:, 1) == weights
 *
 * diagonal and subDiagonal will be modified after the call to this function.
 * diagonal should be N elements
 * subdiagonal should be N-1 elements
 */
blitz::Array<double, 2> CalculateQuadrature(int N, blitz::Array<double, 1> diagonal, blitz::Array<double, 1> subDiagonal);
blitz::Array<double, 2> CalculateQuadrature(int N, PolynomialType type);

} //namespace


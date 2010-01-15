#ifndef ORTHOPOL_TOOLS_H
#define ORTHOPOL_TOOLS_H

#include "../../common.h"

namespace OrthoPol
{

using namespace blitz;

/*
 * Available transformation types
 */
enum PolynomialType
{
	HermitePolynomial = 1,
	LaguerrePolynomial = 2
};

/*
 * Transformation parameters as used by methods in this namespace
 */
struct Parameter
{
	PolynomialType Type;
	double Scaling;
	int HypersphericalRank;
};


/*
 * Computes Associated Laguerre quadrature points and weights by 
 * solving a tridiagonal matrix.
 * 	N  (input): number of gridpoints
 * 	m  (input): the choice of Ass. Laguerre functions, L^(m)
 *	X (output): quadrature points/abscissas
 *	W (output): quadrature weights
 */
void LaguerreQuad(int N, const double &m, Array<double, 1> &X, Array<double, 1> &W);

/*
 * Computes the normalized Associated Laguerre function matrix by
 * recursion.
 *	N  (input): number of gridpoints/functions
 *	m  (input): the choice of Ass. Laguerre functions, L^(m)
 *	X  (input): gridpoints/abscissas
 *	L (output): Ass. Laguerre function matrix
 */
void LaguerreMatrix(int N, const double &m, const Array<double, 1> &X, Array<double, 2> &L);

/*
 * Computes Hermite quadrature points and weights by solving 
 * a tridiagonal matrix.
 * 	N  (input): number of gridpoints
 *	X (output): quadrature points/abscissas
 *	W (output): quadrature weights
 */
void HermiteQuad(int N, Array<double, 1> &X, Array<double, 1> &W);

/*
 * Computes the normalized Hermite function matrix by
 * recursion.
 *	N  (input): number of gridpoints/functions
 *	X  (input): gridpoints/abscissas
 *	H (output): Ass. Laguerre function matrix
 */
void HermiteMatrix(int N, const Array<double, 1> &X, Array<double, 2> &H);

/*
 *	Computes gridpoints and quadrature weights for either
 *	Hermite or Laguerre transforms. 
 * 	N     (input): number of gridpoints
 *	Param (input): transform parametres
 *	X    (output): quadrature points/abscissas
 *	W    (output): quadrature weights
 *
 *  Remarks: The Grid is scaled according to param.Scaling
 */
void ScaledGridAndWeights(int N, const Parameter &param, Array<double, 1> &X, Array<double, 1> &W);
   

/*
 * Calculates all information needed to construct a propagator for the
 * orthogonal polynomial described in param
 */
void OrthoPolSetup(int N, const Parameter &param, Array<double, 1> &eigenvalues, Array<double, 2> &eigenvectors, Array<double, 1> &x, Array<double, 1> &w);

       
}; //Namespace

#endif


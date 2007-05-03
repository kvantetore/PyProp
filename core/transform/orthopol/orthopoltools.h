#ifndef ORTHOPOL_TOOLS_H
#define ORTHOPOL_TOOLS_H

#include "../../common.h"

namespace OrthoPol
{

using namespace blitz;

/*
 * Available transformation types
 */
enum TransformType
{
	HermiteTransform = 1,
	LaguerreTransform = 2
};

/*
 * Transformation parameters as used by methods in this namespace
 */
struct Parameter
{
	TransformType Type;
	double Cutoff;
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
 * Sets up the propagation scheme for the Laguerre transform, 
 * by computing the propagation matrix and the gridpoints.
 *	N      (input): number of gridpoints/functions
 *	cutoff (input): all quadrature points will be less or equal to this
 * 	P     (output): propagation matrix
 */
void SetupLaguerre(int N, const cplx &dt, const double &cutoff, Array<cplx, 2> &P);

/*
 * Sets up the propagation scheme for the Hermite transform, 
 * by computing the propagation matrix and the gridpoints.
 *	N      (input): number of gridpoints/functions
 *	dim    (input): the dimension of the problem
 *	cutoff (input): all quadrature points will be less or equal to this
 * 	P     (output): propagation matrix
 */
void SetupHermite(int N, const cplx &dt, const double &cutoff, Array<cplx, 2> &P);

/*
 *	Computes gridpoints and quadrature weights for either
 *	Hermite or Laguerre transforms. 
 * 	N     (input): number of gridpoints
 *	Param (input): transform parametres
 *	X    (output): quadrature points/abscissas
 *	W    (output): quadrature weights
 */
void GridAndWeights(int N, const Parameter &Param, Array<double, 1> &X, Array<double, 1> &W, double &Alpha);
   
/*
 * Computes Associated Laguerre quadrature points and weights by 
 * solving a tridiagonal matrix.
 * 	N     (input): number of gridpoints
 *	Param (input): transform parametres
 *	dt    (input): time step
 *	P    (output): propagation matrix
 */
void OrthoPolSetup(int N, const Parameter &Param, const cplx &dt, Array<cplx, 2> &P);
       
}; //Namespace

#endif


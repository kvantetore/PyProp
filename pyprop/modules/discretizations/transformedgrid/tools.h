#ifndef TRANSFORMEDGRID_TOOLS_H
#define TRANSFORMEDGRID_TOOLS_H

#include "../../common.h"

namespace TransformedGrid
{
	using namespace blitz;

	/*
	 * Available transformation types
	 */
	enum TransformType
	{
		TransformTypeAlgebraic = 1,      // x(r) = (L - r) / (L + r)
		TransformTypeTrigonometric = 2,  // x(r) = 4/pi * atan(r/L) - 1
		TransformTypeLogarithmic = 3     // x(r) = 1 - 2 * exp( - r / L )
	};

	enum TransformRange
	{
		TransformRangeRadial = 1,        // r = [0, \infty)
		TransformRangeCartesian = 2      // r = (-\infty, \infty)
	};
	
	/*
	 * Transformation parameters as used by methods in this namespace
	 */
	struct Parameter
	{
	    TransformType Type;
		TransformRange Range;
		double Scaling;

		Parameter() : Type(TransformTypeAlgebraic), Range(TransformRangeRadial) {}
	};
	
	/*
	 * Sets up the chebyshev grid and differentiation matrix
	 * The arrays will be resized to the correct size (N+1)
	 * Parameters:
	 * 	  N (input): number of gridpoints
	 * 	  x (output): Chebyshev points
	 * 	  D (output): Chebyshev differentiation matrix
	 */
	void cheb(int N, Array<double, 1> &x, Array<double, 2> &D);

	/*
	 * Sets up the grid and scaling vectors. The arrays will
	 * be resized to the correct size (N)
	 * Parameters:
	 * 	  N (input): Number of gridpoints
	 * 	  param (input): Transformation type
	 *    X (output): X scaling vector (see "the paper" for details)
	 *    Y (output): Y scaling vector (see "the paper" for details)
	 *    r (output): Grid points
	 */
	void XYmat(int N, const Parameter &param, Array<double, 1> &X, Array<double, 1> &Y, Array<double, 1> &r);

	
	/*
	 * Sets up the weights given the radial points. The array weight
	 * will be resized to the correct size (N)
	 * Parameters:
	 * 	  N (input): Number of gridpoints
	 * 	  param (input): Transformation type
	 * 	  r (input): Grid points, as calculated by XYmat
	 * 	  weight (output): integration weights
	 */
	void SetupWeights (int N, const Parameter &param, const Array<double, 1> &r, Array<double, 1> &weight);

	/*
	 * Sets up a propagation scheme, with forward and inverse transforms, as well as 
	 * the eigenvalues of the propagation operator
	 * Parameters:
	 *    N (input): Number of gridpoints
	 *    transformType (input): which transformation to choose
	 *    	    1: x(r) = (L - r) / (L + r)
	 *    	    2: x(r) = 4/pi * atan(r/L) - 1
	 *    	    3: x(r) = 1 - 2 * exp( - r / L )
	 *   		(tore: i believe the parameter L is fixed to 1 for
	 *   		all transforms at the moment)
	 *    WR (output): Eigenvalues of the operator
	 *    VR (output): Eigenvectors of the operator
	 *    INVV (output): Inverse eigenvectors of the operator
	 */
	void setup (int N, const Parameter &param, Array<double, 1> &WR, Array<double,2> &VR, Array<double,2> &INVV);
		
};

#endif


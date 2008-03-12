#include "bspline.h"
#include "../../utility/blitzblas.h"

namespace BSpline
{

/*! \fn double EvaluateBSpline(double x, int splineOrder, int leftKnotPoint)
 * Evaluate B-spline value B(x) at point 'x' over knot sequence. B-spline
 * is of order 'splineOrder'. B(x) is non-zero over the interval defined
 * by 'leftKnotPoint' and 'leftKnotPoint' + 'splineOrder'.
 * 
 * The calculation is performed using the three-point recursion formula
 * for splines,
 *		
 *     B(k, i, x) = (x - t_i) / (t_(i+k-1) - t_i) * B(k-1, i, x)
 *                  + (t_(i+k) - x) / (t_(i+k) - t(i+1)) * B(k-1, i+1, x)
 *
 * where it is understood that 
 *
 *     B(1, i, x) = 1  , t_i <= x < t_(i+1)
 *     B(1, i, x) = 0  , otherwise
 * 
 */
double BSpline::EvaluateBSpline(double x, int splineOrder, int leftKnotPoint)
{
	
	// Return value
	double B1_0 = 0.0;
		
	if (splineOrder == 1)
	{ 

		/*
		 * If required spline is of order 1, check requirement on knot points and x,
		 * returning B = 1 if all is good, otherwise 0.
		 */
		if ((KnotSequence(leftKnotPoint) <= x) && (x < KnotSequence(leftKnotPoint + 1)))
		{
			B1_0 = 1.0;
		} 
		else 
		{
			B1_0 = 0.0;
		}
	}
	else
	{
		/*
		 * Compute b-spline of order > 1. We must check that b-spline values returned 
		 * from next recursion call is not zero, in which case we avoid 0/0-problem
		 * in recursion formula by removing offending term.
		 */
		double startKnot_0 = KnotSequence(leftKnotPoint);
		double endKnot_0 = KnotSequence(leftKnotPoint + splineOrder - 1);
		double startKnot_1 = KnotSequence(leftKnotPoint + 1);
		double endKnot_1 = KnotSequence(leftKnotPoint + splineOrder);

		//Recursion
		double B0_0 = EvaluateBSpline(x, splineOrder - 1, leftKnotPoint);
		double B0_1 = EvaluateBSpline(x, splineOrder - 1, leftKnotPoint + 1);

		if (B0_0 > eps)
		{
			B1_0 += (x - startKnot_0) / (endKnot_0 - startKnot_0) * B0_0;
		}
		if (B0_1 > eps)
		{
			B1_0 += (endKnot_1 - x) / (endKnot_1 - startKnot_1) * B0_1;
		}
	}

	return B1_0;

} // End function EvaluateBSpline


/*
 * Evaluate second derivate of B-spline using recursive formula 
 */
double BSpline::EvaluateBSplineDerivative2(double x, int splineOrder, int leftKnotPoint)
{
	
	// Return value
	double bspline = 0.0;

	// Shorthands
	int t = leftKnotPoint;
	int k = splineOrder;

	double B1 = EvaluateBSpline(x, splineOrder - 2, leftKnotPoint);
	double B2 = EvaluateBSpline(x, splineOrder - 2, leftKnotPoint + 1);
	double B3 = EvaluateBSpline(x, splineOrder - 2, leftKnotPoint + 2);

	if (B1 > eps)
	{
		double b_ = (k - 2.0) / ( KnotSequence(t + k - 2) - KnotSequence(t) );
		B1 *= b_;
	}
	
	if (B2 > eps)
	{
		double c_ = (k - 2.0) / ( KnotSequence(t + k - 1) - KnotSequence(t + 1));
		B2 *= c_;
	}

	if (B3 > eps)
	{
		double e_ = (k  - 2.0) / ( KnotSequence(t + k) - KnotSequence(t + 2) );
		B3 *= e_;
	}

	if ((B1 > eps) || (B2 > eps))
	{
		double a_ = (k - 1.0) / ( KnotSequence(t + k - 1) - KnotSequence(t) );
		bspline += a_ * (B1 - B2);
	}

	if ((B2 > eps) || (B3 > eps))
	{
		double d_ = (k - 1.0) / ( KnotSequence(t + k) - KnotSequence(t + 1) );
		bspline -= d_ * (B2 - B3);
	}

	return bspline;

} // End function EvaluateBSplineDerivative2


/*
 * Tabulate b-spline values at quadrature points. This is done
 * to speed up later computations, where we avoid repeated calls
 * to the generating recursion function.
 */
void BSpline::CreateBSplineTable()
{
	int numberOfQuadPoints = Nodes.extent(0);
	int xSize = numberOfQuadPoints * MaxSplineOrder;
	
	// Rescale arrays
	QuadratureGrid.resize(NumberOfBSplines, xSize);
	BSplineTable.resize(NumberOfBSplines, xSize);
	ScaledWeights.resize(NumberOfBSplines, xSize);
	Ones.resize(xSize);

	QuadratureGrid = 0.0;
	BSplineTable = 0.0;
	ScaledWeights = 0.0;
	Ones = 1.0;

	for (int i = 0; i < NumberOfBSplines; i++)
	{
		int xIndex = 0;

		for (int j = 0; j < MaxSplineOrder; j++)
		{
			double a = KnotSequence(i+j);
			double b = KnotSequence(i+j+1);
			if( fabs(a - b) < eps) { continue; }

			for (int k = 0; k < numberOfQuadPoints; k++)
			{
				double x = ScaleAndTranslate(Nodes(k), a, b);
				QuadratureGrid(i, xIndex) = x;
				ScaledWeights(i, xIndex) = (b - a) / 2.0 * Weights(k);
				BSplineTable(i, xIndex) = EvaluateBSpline(x, MaxSplineOrder, i);
				xIndex += 1;
			}
		}
	}
	
} // End function CreateBSplineTable()


blitz::Array<double, 1> BSpline::GetBSpline(int bsplineIndex)
{
	//int xSize = Nodes.extent(0) * MaxSplineOrder;
	return BSplineTable( bsplineIndex, blitz::Range::all() );
}


blitz::Array<double, 1> BSpline::GetBSplineDerivative2(int bsplineIndex)
{
	using namespace blitz;
	VectorType D2B = BSplineDerivative2Table( bsplineIndex, Range::all() );
	return D2B;
}


/*
 * Tabulate b-spline second derivative values at quadrature points. 
 * This is done to speed up later computations, where we avoid 
 * repeated calls to the generating recursion function.
 */
void BSpline::CreateBSplineDerivative2Table()
{
	int numberOfQuadPoints = Nodes.extent(0);
	int xSize = numberOfQuadPoints * MaxSplineOrder;
	
	// Rescale arrays
	BSplineDerivative2Table.resize(NumberOfBSplines, xSize);

	BSplineDerivative2Table = 0.0;

	for (int i = 0; i < NumberOfBSplines; i++)
	{
		int xIndex = 0;

		for (int j = 0; j < MaxSplineOrder; j++)
		{
			double a = KnotSequence(i+j);
			double b = KnotSequence(i+j+1);
			if( fabs(a - b) < eps) { continue; }

			for (int k = 0; k < numberOfQuadPoints; k++)
			{
				double x = ScaleAndTranslate(Nodes(k), a, b);
				BSplineDerivative2Table(i, xIndex) = 
					EvaluateBSplineDerivative2(x, MaxSplineOrder, i);
				xIndex += 1;
			}
		}
	}
	
} // End function CreateBSplineDerivative2Table()


/*
 * Same as CreateBSplineTable, but storage in
 * BLAS order.
 */ 
/*void BSpline::CreateBSplineTableBlas()
{
	int numberOfQuadPoints = QuadratureGridGlobal.extent(0);
	//int xSize = numberOfQuadPoints * MaxSplineOrder;
	int k = MaxSplineOrder;
	int numberOfNodes = Nodes.extent(0);
	int numberOfBands = (2 * k - 1) * numberOfNodes;
	
	// Rescale arrays
	BSplineTableBlas.resize(numberOfQuadPoints, numberOfBands);

	BSplineTableBlas = 0.0;

	for (int j = 0; j < numberOfQuadPoints; j++)
	{
		int iMin = std::max(j - numberOfBands + 1, 0);
		int iMax = std::min(j + numberOfBands, NumberOfBSplines);
		for (int i = iMin; i < iMax; i++)
		{
			int I = k - j + i;
			int J = j;
			double x = QuadratureGridGlobal(j);
			//BSplineTableBlas(J,I) = BSplineTable(i, xIndex)
			BSplineTableBlas(J,I) = EvaluateBSpline(x, k, i);
		}
	}
	
} // End function CreateBSplineTableBlas	
*/


/*
 * Compute overlap integral between b-splines B_i and B_j.
 * Due to the compact support of b-splines, a given b-spline
 * will only overlap with k other splines (k: spline order).
 * This routine accounts for this and returns zero if |i-j| >= k.
 * NOTE: For now we must have j >= i!
 */
double BSpline::BSplineOverlapIntegral(VectorType f, int i, int j)
{
	using namespace blitz;

	// Overlap only nonzero if |i - j| < splineOrder
	if ( abs(i - j) >= MaxSplineOrder )
	{
		return 0.0;
	}
	else
	{
		//
		// Since b-splines are only evaluated and stored on the part
		// of the grid where they are non-zero, we need to figure out
		// how much of the two given b-splines overlap (i.e. indices)
		// and only sum over this area.
		//
		int leftIndex = i;
		int rightIndex = j;
		int startIndex, stopIndex;
		if ( fabs(KnotSequence(i) - KnotSequence(j)) < eps )
		{
			startIndex = 0;
			stopIndex = BSplineTable.extent(1) - 1;
		}
		else
		{
			startIndex = ComputeStartIndex(leftIndex, rightIndex);
			stopIndex = ComputeStopIndex(startIndex, rightIndex);
		}

		// Slice left and right b-splines + function and scaled integration weights
		VectorType B_left = BSplineTable(leftIndex, Range(startIndex, toEnd));
		VectorType B_right = BSplineTable(rightIndex, Range(fromStart, stopIndex));
		VectorType ScaledWeightsSlice = ScaledWeights(j, Range(0, stopIndex));
		VectorType fSlice = f(Range(startIndex, toEnd));

		// Calculate overlap
		double I = 0;
		for (int i = 0; i < fSlice.extent(0); i++)
		{
			I += ScaledWeightsSlice(i) * fSlice(i) * B_left(i) * B_right(i);
		}
		//double I = sum(ScaledWeightsSlice * fSlice * B_left * B_right);

		return I;
	}

	return 0.0;
}


/*
 * Compute overlap integral between b-splines B_i d2/dx2 B_j.
 * Due to the compact support of b-splines, a given b-spline
 * will only overlap with 2k-1 other splines (k: spline order).
 * This routine accounts for this and returns zero if |i-j| >= k.
 * NOTE: For now we must have j >= i!
 */
double BSpline::BSplineDerivative2OverlapIntegral(int i, int j)
{
	using namespace blitz;

	// Overlap only nonzero if |i - j| < splineOrder
	if ( abs(i - j) >= MaxSplineOrder )
	{
		return 0.0;
	}
	else
	{
		//
		// Since b-splines are only evaluated and stored on the part
		// of the grid where they are non-zero, we need to figure out
		// how much of the two given b-splines overlap (i.e. indices)
		// and only sum over this area.
		//
		int leftIndex = i;
		int rightIndex = j;
		int startIndex, stopIndex;
		if ( fabs(KnotSequence(i) - KnotSequence(j)) < eps )
		{
			startIndex = 0;
			stopIndex = BSplineTable.extent(1) - 1;
		}
		else
		{
			startIndex = ComputeStartIndex(leftIndex, rightIndex);
			stopIndex = ComputeStopIndex(startIndex, rightIndex);
		}

		// Slice left and right b-splines + function and scaled integration weights
		VectorType B_left = BSplineTable(leftIndex, Range(startIndex, toEnd));
		VectorType B_right = BSplineDerivative2Table(rightIndex, Range(fromStart, stopIndex));
		VectorType ScaledWeightsSlice = ScaledWeights(j, Range(0, stopIndex));

		// Calculate overlap
		double I = sum(ScaledWeightsSlice * B_left * B_right);

		return I;
	}

	return 0.0;
}

cplx BSpline::ProjectOnBSpline(VectorTypeCplx func, int i)
{
	using namespace blitz;
	return sum( ScaledWeights(i, Range::all()) * func * BSplineTable( i, Range::all() ) );
}


/*
 * Evaluate a B-spline on given grid points
 */
blitz::Array<double, 1> BSpline::EvaluateBSplineOnGrid(VectorType grid, int bSplineNumber)
{
	int gridSize = grid.extent(0);
	VectorType bspline = VectorType(gridSize);
	for (int i = 0; i < gridSize; i++)
	{
		double x = grid(i);
		bspline(i) = EvaluateBSpline(x, MaxSplineOrder, bSplineNumber);
	}

	return bspline;
}


/*! \fn blitz::Array<cplx, 1> ExpandFunctionInBSplines(object func)
 * Expand a function f in the B-spline basis. Since B-splines are not
 * orthogonal, the usual projection method yields a linear system of
 * equations, Sc = b, where S is the B-spline overlap matrix, c are the 
 * expansion coefficients and b is the projection of f onto the basis.
 * The matrix S is banded with 2k - 1 bands (k is B-spline order).
 * Inverting S yields the expansion coefficients, c = S'b.
 */
blitz::Array<cplx, 1> BSpline::ExpandFunctionInBSplines(object func)
{
	using namespace blitz;

	VectorTypeCplx b = VectorTypeCplx(NumberOfBSplines);
	int gridChunkSize = QuadratureGrid.extent(1);
	VectorTypeCplx functionGridValues = VectorTypeCplx(gridChunkSize);

	for (int i = 0; i < NumberOfBSplines; i++)
	{
		for (int j = 0; j < gridChunkSize; j++)
		{
			functionGridValues(j) = extract<cplx>(func(QuadratureGrid(i, j)));
		}

		b(i) = ProjectOnBSpline(functionGridValues, i);
	}

	/* 
	 * Solving banded linear system of equations 
	 *  to obtain expansion coefficients 
	 */
	SetupOverlapMatrixFull();
	int offDiagonalBands = MaxSplineOrder - 1;
	VectorTypeInt pivot(NumberOfBSplines);
	pivot = 0;
	
	linalg::LAPACK<cplx> lapack;
	lapack.SolveGeneralBandedSystemOfEquations(OverlapMatrixFull, pivot, b, 
		offDiagonalBands, offDiagonalBands);

	OverlapMatrixFullComputed = false;
	
	return b;
}

/* 
 * Expand a function f in the B-spline basis. The function is precalculated
 * on the quadrature grid, and passed as a complex 1D blitz::Array.
 *
 * There are two algorithms implemented here:
 *
 *     0. Sum is performed with explicit loop. Performance: approx. two
 *        times faster than 1.
 *     1. Sum is performed with blitz++ sum-function (which is slow).
 */
void BSpline::ExpandFunctionInBSplines(blitz::Array<cplx, 1> input, blitz::Array<cplx, 1> output)
{
	using namespace blitz;

	int gridChunkSize = QuadratureGrid.extent(1);
	int globalGridSize = QuadratureGridGlobal.extent(0);

	if (ProjectionAlgorithm == 0)
	{
		for (int i = 0; i < NumberOfBSplines; i++)
		{
			int startIndex = GetGridIndex(i);
			int stopIndex = std::min(startIndex + gridChunkSize - 1, globalGridSize - 1);
			int curChunkSize = stopIndex - startIndex;

			cplx tempValue = 0;
			for (int j=0; j<curChunkSize; j++)
			{
				tempValue += input(j+startIndex) * ScaledWeights(i, j) * BSplineTable(i, j);
			}
			output(i) = tempValue;

		}
	}

	else if (ProjectionAlgorithm == 1)
	{
		VectorTypeCplx functionSlice = VectorTypeCplx(gridChunkSize);
		for (int i = 0; i < NumberOfBSplines; i++)
		{
			int startIndex = GetGridIndex(i);
			int stopIndex = std::min(startIndex + gridChunkSize - 1, globalGridSize - 1);
			functionSlice = 0.0;
			functionSlice( Range(0, stopIndex-startIndex) ) = 
				input( Range(startIndex, stopIndex) );

			// Compute projection on b-spline i
			output(i) = ProjectOnBSpline(functionSlice, i);
		}
	}

	else

	{
		cout << "ERROR: Unknown projection algorithm " << ProjectionAlgorithm << endl;
	}
 	
/*	else if (ProjectionAlgorithm == 2)
	{
		double* weightData = &ScaledWeights(0,0);
		double* splineData = &BSplineTable(0,0);
		cplx* outputData = &output(0);
		int inputStride = input.stride(0);

		for (int i = 0; i < NumberOfBSplines; i++)
		{
			int startIndex = GetGridIndex(i);
			int stopIndex = std::min(startIndex + gridChunkSize - 1, globalGridSize - 1);
			int curChunkSize = stopIndex - startIndex;
			
			cplx* inputData = & input(startIndex);


			for (int j=0; j<curChunkSize; j++)
			{
				(*outputData) += (*inputData) * (*weightData) * (*splineData);
				inputData += inputStride;
				weightData++;
				splineData++;
			}

			outputData++;
		} 

	}
*/

	/* 
	 * Solving banded linear system of equations 
	 *  to obtain expansion coefficients 
	 */
	SetupOverlapMatrixFull();
	int offDiagonalBands = MaxSplineOrder - 1;
	VectorTypeInt pivot(NumberOfBSplines);
	pivot = 0;
	
	linalg::LAPACK<cplx> lapack;
	lapack.SolveGeneralBandedSystemOfEquations(OverlapMatrixFull, pivot, output, 
		offDiagonalBands, offDiagonalBands);

	OverlapMatrixFullComputed = false;
}


/*
 * From a given sequence of b-spline coefficients, construct a function on given grid.
 */
blitz::Array<cplx, 1> BSpline::ConstructFunctionFromBSplineExpansion(VectorTypeCplx c, VectorType grid)
{
	int gridSize = grid.extent(0);
	int numberOfCoeffs = c.extent(0);
	VectorTypeCplx f = VectorTypeCplx(gridSize);
	f = 0.0;

	for (int i = 0; i < numberOfCoeffs; i++)
	{
		f += c(i) * EvaluateBSplineOnGrid(grid, i);	
	}
	
	return f;
}

/*
 * Reconstruct a function on quadrature grid points
 */
blitz::Array<cplx, 1> BSpline::ConstructFunctionFromBSplineExpansion(VectorTypeCplx c)
{
	using namespace blitz;

	// Size of global quadrature grid
	int gridSize = QuadratureGridGlobal.extent(0);

	// Size of grid slice where a b-spline is defined
	int gridSliceSize = QuadratureGrid.extent(1);

	int breakpointSize = BreakpointSequence.extent(0);

	// This will hold the reconstructed function
	Array<cplx, 1> function(gridSize);
	function = 0.0;

	// Index on actual (global) quadrature grid
	int gridIndex = 0;

	// B-spline index
	int bSplineIndex = 0;

	for (int i = 0; i < breakpointSize - 1; i++)
	{
		
		int knotDegeneracy = MaxSplineOrder - ContinuitySequence(i);
		for (int j = 0; j < knotDegeneracy; j++)
		{
			/*
			 * Need to figure out where the grid slice indices on which a given
			 * b-spline is evaluated. Degenerate knot points are a particular pain.
			 * However, since all b-splines are tabulated on an equal number of point
			 * (storing zeros outside their nonzero region), we can ignore the fact
			 * that they only are nonzero on less than MaxSplineOrder breakpoints.
			 * StopIndex for grid slice is therefore just startIndex + gridSliceSize,
			 * except at the end of the grid, where we stop at (global) gridSize.
			 * Likewise the k last b-splines must be sliced shorter to match the
			 * (shorter) length of the grid slice.
			 */
			int startIndex = gridIndex;
			int stopIndex = std::min(startIndex + gridSliceSize, gridSize);
			int bSplineGridStopIndex = std::min(gridSliceSize, stopIndex - startIndex);

			function( Range(startIndex, stopIndex - 1) ) += c(bSplineIndex)
				* BSplineTable(bSplineIndex, Range(fromStart, bSplineGridStopIndex - 1));
			
			bSplineIndex++;
			if (bSplineIndex >= NumberOfBSplines) { return function; }
		}

		gridIndex += Nodes.extent(0);
	}

	return function;
}


/*
 * Compute b-spline overlap matrix. The matrix is banded due to
 * compactness of the b-splines, and symmetric positive definite.
 * For easier use with LAPACK routines, we store only the upper
 * bands in the manner specified by LAPACK (f.ex. zpbsvx).
 *
 * LAPACK indices are computed from
 *
 *     I =  MaxSplineOrder - 1 + i - j
 *     J = j
 *	
 *	where we have -1 from FORTRAN->C conversion of original
 *	formula.
 */
void BSpline::ComputeOverlapMatrix()
{
	OverlapMatrix.resize(NumberOfBSplines, MaxSplineOrder);
	OverlapMatrix = 0.0;

	for (int i = 0; i < NumberOfBSplines; i++)
	{
		int jMax = std::min(i+MaxSplineOrder, NumberOfBSplines);
		for (int j = i; j < jMax ; j++)
		{
			OverlapMatrix(j, MaxSplineOrder - 1 + i - j) = BSplineOverlapIntegral(i, j);
		}
	}

	OverlapMatrixComputed = true;
}


/*
 * Setup "full" overlap matrix for use with LAPACK routine zgbsv.
 * Both upper and lower bands are stored, in a manner required by
 * said routine.
 */
void BSpline::SetupOverlapMatrixFull()
{
	int k = MaxSplineOrder;

	if (!OverlapMatrixComputed)
	{ 
		ComputeOverlapMatrix();
	}
	
	if (OverlapMatrixFullComputed) { return; }

	if (OverlapMatrixFull.extent(0) != NumberOfBSplines)
	{
		int ldab = 3 * k - 2;
		OverlapMatrixFull.resize(NumberOfBSplines, ldab);
	}

	for (int i = 0; i < NumberOfBSplines; i++)
	{
		int jMax = std::min(i+MaxSplineOrder, NumberOfBSplines);
		for (int j = i; j < jMax ; j++)
		{
			// OverlapMatrix is stored according to (another) LAPACK
			// specification, so first we must map its indices to
			// "ordinary" matrix indices.
			int J = j;
			int I = k - 1 + i - j;

			// Then we use the zgbsv-map. However, since this is
			// a full matrix we must transpose to obtain lower
			// triangle. This is done by flipping (i,j)->(j,i)
			// where (i,j) are the "ordinary" indices.
			int J2 = j;
			int I2 = 2 * k - 2 + i - j; 
			int J_lower = i;
			int I_lower = 2 * k - 2 + j - i; 
			OverlapMatrixFull(J2, I2) = OverlapMatrix(J, I);
			OverlapMatrixFull(J_lower, I_lower) = OverlapMatrix(J, I);
		}
	}

	OverlapMatrixFullComputed = true;
}


/*
 * Setup overlap matrix for use with BLAS routine zgbsv.
 * Both upper and lower bands are stored, in a manner required by
 * said routine.
 */
void BSpline::SetupOverlapMatrixBlas()
{
	int k = MaxSplineOrder;
	int N = NumberOfBSplines;

	//Check that overlap matrix has been computed
	if (!OverlapMatrixComputed)
	{ 
		ComputeOverlapMatrix();
	}

	//Check if BLAS overlap matrix has been set up
	if (OverlapMatrixBlasComputed) { return; }

	//Check that size of OverlapMatrixBlas is correct; otherwise resize
	if ( (OverlapMatrixBlas.extent(0) != N) || (OverlapMatrixBlas.extent(1) != k) )
	{
		OverlapMatrixBlas.resize(N, k);
	}

	
	OverlapMatrixBlas = 0;

	for (int j = 0; j < N; j++)
	{
		int iMax = std::min(j + k, N);
		for (int i = j; i < iMax; i++)
		//for (int i = j; i < N; i++)
		{
			/* 
			 * OverlapMatrix is stored according to (another) LAPACK
			 * specification, so first we must map its indices to
			 * "ordinary" matrix indices.
			 */
			int J = i;
			int I = k - 1 + j - i;

			//BLAS index map from "normal" indices
			int Jb = j;
			int Ib = i - j;

			//Set up blas matrix
			OverlapMatrixBlas(Jb, Ib) = OverlapMatrix(J, I);
			//OverlapMatrixBlas(Jb, Ib) = BSplineOverlapIntegral(j, i);
		}
	}

	OverlapMatrixBlasComputed = true;
}


} //Namespace


#ifndef BSPLINES_H
#define BSPLINES_H

#include "../../utility/blitzlapack.h"
#include "../../utility/boostpythonhack.h"

namespace BSpline
{

class BSpline
{

public:
	typedef shared_ptr<BSpline> Ptr;
	typedef blitz::Array<double, 1> VectorType;
	typedef blitz::Array<cplx, 1> VectorTypeCplx;
	typedef blitz::Array<double, 2> MatrixType;
	typedef blitz::Array<cplx, 2> MatrixTypeCplx;
	typedef blitz::Array<int, 1> VectorTypeInt;

private:

	double eps;
	
	VectorType BreakpointSequence;
	VectorTypeInt ContinuitySequence;
	VectorType KnotSequence;
	VectorTypeInt KnotGridIndexMap;    // Map a knot index to corresponding global grid index
	VectorTypeInt TopKnotMap;          // Map a knot index to "top" knot index
	VectorType Weights;
	VectorType Nodes;
	VectorType Ones;
	VectorType QuadratureGridGlobal;    // Repeated concatenation of quadrature points
	VectorType QuadratureWeightsGlobal; // Repeated concatenation of quadrature weights

	MatrixType BSplineTable;
	MatrixTypeCplx BSplineTableBlas;
	MatrixType BSplineDerivative2Table;
	MatrixType QuadratureGrid;
	MatrixType ScaledWeights;
	MatrixType OverlapMatrix;
	MatrixTypeCplx OverlapMatrixFull;
	MatrixType OverlapMatrixBlas;

	//Needed by zpbsvx expert interface
	MatrixTypeCplx OverlapMatrixExpert;
	MatrixTypeCplx FactoredOverlapMatrix;
	VectorType OverlapScaleFactors;
	VectorType ForwardErrorEstimate;
	VectorType BackwardErrorEstimate;

	bool OverlapMatrixComputed;
	bool OverlapMatrixFullComputed;
	bool OverlapMatrixBlasComputed;

public:

	// Data members	
	int NumberOfBSplines;
	int MaxSplineOrder;

	int ProjectionAlgorithm;
	int LapackAlgorithm;

	// Constructor
	BSpline() 
	{ 
		eps = 1e-15; 
		OverlapMatrixComputed = false;
		OverlapMatrixFullComputed = false;
		OverlapMatrixBlasComputed = false;
		ProjectionAlgorithm = 0;
		LapackAlgorithm = 0;
	}

	// Destructor
	virtual ~BSpline() {}

	// Functions to get grid-related sequences
	VectorType GetBreakpointSequence() { return BreakpointSequence; }
	VectorTypeInt GetContinuitySequence() { return ContinuitySequence; }
	VectorType GetKnotSequence() { return KnotSequence; }
	VectorTypeInt GetKnotGridIndexMap() { return KnotGridIndexMap; }
	VectorTypeInt GetTopKnotMap() { return TopKnotMap; }
	VectorType GetWeights() { return Weights; }
	VectorType GetNodes() { return Nodes; }
	VectorType GetQuadratureGrid(int i) { return QuadratureGrid(i, blitz::Range::all()); }
	MatrixType GetBSplineOverlapMatrix() { return OverlapMatrix; }
	MatrixTypeCplx GetBSplineOverlapMatrixFull() { SetupOverlapMatrixFull(); return OverlapMatrixFull; }
	MatrixType GetBSplineOverlapMatrixBlas() { SetupOverlapMatrixBlas(); return OverlapMatrixBlas; }
	MatrixTypeCplx GetBSplineOverlapMatrixExpert() { return OverlapMatrixExpert; }
	VectorType GetQuadratureGridGlobal() { return QuadratureGridGlobal; }
	VectorType GetQuadratureWeightsGlobal() { return QuadratureWeightsGlobal; }

	// Functions to set the length of various sequences
	void ResizeBreakpointSequence(int size) { BreakpointSequence.resize(size); }
	void ResizeContinuitySequence(int size) { ContinuitySequence.resize(size); }
	void ResizeKnotSequence(int size) { KnotSequence.resize(size); }
	void ResizeKnotGridIndexMap(int size) { KnotGridIndexMap.resize(size); }
	void ResizeTopKnotMap(int size) { TopKnotMap.resize(size); }
	void ResizeWeights(int size) { Weights.resize(size); }
	void ResizeNodes(int size) { Nodes.resize(size); }
	void ResizeQuadratureGridGlobal(int size) { QuadratureGridGlobal.resize(size); }
	void ResizeQuadratureWeightsGlobal(int size) { QuadratureWeightsGlobal.resize(size); }
	//void ResizeOverlapMatrix(int size) { OverlapMatrix.Resize(size); }

	// B-spline evaluation functions
	double EvaluateBSpline(double, int, int);
	double EvaluateBSplineDerivative2(double, int, int);
	VectorType EvaluateBSplineOnGrid(VectorType, int);
	VectorType GetBSpline(int);
	VectorType GetBSplineDerivative2(int);
	double BSplineOverlapIntegral(int i, int j) { return BSplineOverlapIntegral(Ones, i, j); }
	double BSplineOverlapIntegral(VectorType, int, int);
	cplx ProjectOnBSpline(VectorTypeCplx, int); 
	//double EvaluateBSplineDerivative(int int);
	double BSplineDerivative2OverlapIntegral(int, int);
	
	void SolveForOverlapMatrix(VectorTypeCplx vector);

	template<class TBase> TBase BSplineGlobalOverlapIntegral(blitz::Array<TBase, 1> func, int derivative, int i, int j);
	//For python wrapping
	double BSplineGlobalOverlapIntegral_double(blitz::Array<double, 1> func, int derivative, int i, int j)
	{
		return BSplineGlobalOverlapIntegral(func, derivative, i, j);
	}

	cplx BSplineGlobalOverlapIntegral_cplx(blitz::Array<cplx, 1> func, int derivative, int i, int j)
	{
		return BSplineGlobalOverlapIntegral(func, derivative, i, j);
	}

	void CreateBSplineTable();
	//void CreateBSplineTableBlas();
	void CreateBSplineDerivative2Table();
	void ComputeOverlapMatrix();
	void SetupOverlapMatrixFull();
	void SetupOverlapMatrixBlas();
	void SetupOverlapMatrixExpert();

	// B-spline-expansion related functions
	VectorTypeCplx ExpandFunctionInBSplines(object);
	
	// Expand pre-evaluated function in b-spline basis
	void ExpandFunctionInBSplines(blitz::Array<cplx, 1> input, blitz::Array<cplx, 1> output); 	

	// Reconstruct on quadrature points
	VectorTypeCplx ConstructFunctionFromBSplineExpansion(VectorTypeCplx);

	// Reconstruct on given grid
	VectorTypeCplx ConstructFunctionFromBSplineExpansion(VectorTypeCplx, VectorType); 
	
	// Various small convenient functions
	int ComputeStartIndex(int, int);
	int ComputeStopIndex(int, int);
	int ComputeIndexShift(int);
	int GetGridIndex(int);
	double ScaleAndTranslate(double, double, double);
	
};


/*
 * Inline functions
 */
inline int BSpline::ComputeStartIndex(int leftIndex, int rightIndex)
{
	int shift = ComputeIndexShift(leftIndex);
	int xSize = Nodes.extent(0);
	return xSize * (rightIndex - leftIndex - shift);
}


inline int BSpline::ComputeStopIndex(int startIndex, int rightIndex)
{
	int xSize = Nodes.extent(0);
	int shift = ComputeIndexShift(rightIndex);
	return BSplineTable.extent(1) - startIndex - xSize * shift - 1;
}


inline int BSpline::ComputeIndexShift(int i)
{
	return TopKnotMap(i) - i;
}

/*
 * Scale and translate a coordinate from [-1,1] to [a,b]
 */
inline double BSpline::ScaleAndTranslate(double x, double a, double b)
{
	return (x - 1.0) * (b - a) / 2.0 + b;
}


/*
 * Get the global grid point correspoding to a knot point
 */
inline int BSpline::GetGridIndex(int knotIndex)
{
	return KnotGridIndexMap(knotIndex);
}



}; //Namespace

#endif


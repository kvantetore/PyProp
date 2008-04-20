#ifndef OVERLAPMATRIX_H
#define OVERLAPMATRIX_H

#include "../common.h"

class OverlapMatrixEvaluator
{
public: 
	OverlapMatrixEvaluator() {}
	virtual ~OverlapMatrixEvaluator() {}

	virtual double operator()(int row, int col) const = 0;
};

class OverlapMatrixEvaluatorDiagonal : public OverlapMatrixEvaluator
{
public: 
	OverlapMatrixEvaluatorDiagonal(blitz::Array<double, 1> weights) : Weights(weights) {}
	virtual ~OverlapMatrixEvaluatorDiagonal() {}

	virtual double operator()(int row, int col) const 
	{
		if (row != col) 
		{
			return 0;
		}
		else 
		{
			return Weights(row);
		}
	}

private:
	blitz::Array<double, 1> Weights;
};


/*
 * Class to deal with overlap matrices. 
 *
 * BasisSize is the number of elements in the basis (NumberOfBSplines)
 * SuperDiagonals is the number of superdiagonals in the overlap matrix (MaxSplineOrder - 1)
 *
 */
class OverlapMatrix
{
public:
	typedef boost::shared_ptr<OverlapMatrix> Ptr;
	typedef blitz::Array<cplx, 2> MatrixType;
	typedef blitz::Array<double, 2> MatrixTypeReal;
	typedef blitz::Array<cplx, 1> VectorType;
	typedef blitz::Array<cplx, 3> TensorType;

	OverlapMatrix(int basisSize, int superDiagonals, const OverlapMatrixEvaluator &evaluator) :
		SuperDiagonals(superDiagonals),
		BasisSize(basisSize)
	{ 
		Setup(evaluator);
	}

	MatrixType GetOverlapHermitianUpper()
	{
		return OverlapHermitianUpper;
	}

	MatrixType GetOverlapHermitianLower()
	{
		return OverlapHermitianLower;
	}

	MatrixType GetOverlapCholeskyUpper()
	{
		return OverlapCholeskyUpper;
	}

	MatrixTypeReal GetOverlapFullRow()
	{
		return OverlapFullRow;
	}

	MatrixTypeReal GetOverlapFullCol()
	{
		return OverlapFullCol;
	}

	int GetSuperDiagonals()
	{
		return SuperDiagonals;
	}

	int GetBasisSize()
	{
		return BasisSize;
	}

	void Setup(const OverlapMatrixEvaluator &evaluator);
	//Vector
	void MultiplyOverlapVector(const VectorType &source, VectorType &dest);
	void MultiplyOverlapVector(VectorType &vector);
	void SolveOverlapVector(VectorType &vector);
	//Tensor is a 3D Array, with the rank of iterest as rank 1
	void MultiplyOverlapTensor(const TensorType &source, TensorType &dest);
	void MultiplyOverlapTensor(TensorType &vector);
	void SolveOverlapTensor(TensorType &vector);

	/*
	 * Sqrt-routines apply only the square root of the overlap matrix S.
	 * as we have the Cholesky factorization S = R* R, this amounts to only applying 
	 * one of the Cholesky factors
	 *
	 * *Conj routines apply conj(R)
	 */
	//Vector
	void MultiplySqrtOverlapVector(bool conjugate, VectorType &vector);
	void SolveSqrtOverlapVector(bool conjugate, VectorType &vector);
	//Tensor
	void MultiplySqrtOverlapTensor(bool conjugate, TensorType &vector);
	void SolveSqrtOverlapTensor(bool conjugate, TensorType &vector);

private:
	int SuperDiagonals;
	int BasisSize;

	MatrixType OverlapHermitianUpper;
	MatrixType OverlapHermitianLower;
	MatrixType OverlapCholeskyUpper;

	//Used by inner product
	MatrixTypeReal OverlapFullRow;
	MatrixTypeReal OverlapFullCol;
};

#endif


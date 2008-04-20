#include "overlapmatrix.h"

#include "../utility/blitzlapack.h"
#include "../utility/blitztricks.h"
#include "../krylov/piram/piram/blitzblas.h"

using namespace blitz::linalg;

void OverlapMatrix::Setup(const OverlapMatrixEvaluator &evaluator)
{
	/*  From BLAS documentation
	 *  upper hermitian banded storage
	 *
	 *  DO J = 1, N
	 *      M = K + 1 - J
	 *      DO I = MAX( 1, J - K ), J
	 *          A( M + I, J ) = matrix( I, J )
	 *      ENDDO
	 *  ENDDO
	 */

	//Hermitian Upper storage
	OverlapHermitianUpper.resize(BasisSize, SuperDiagonals+1);
	OverlapHermitianUpper = 0.0;
	for (int i=0; i<BasisSize; i++)
	{
		int jMax = std::min(i+SuperDiagonals+1, BasisSize);
		for (int j = i; j < jMax ; j++)
		{
			OverlapHermitianUpper(j, SuperDiagonals + i - j) = evaluator(i, j);
		}
	}

	//Factorize the overlap matrix
	OverlapCholeskyUpper.reference( OverlapHermitianUpper.copy() );
	blitz::linalg::LAPACK<cplx> lapack;
	lapack.CalculateCholeskyFactorizationPositiveDefiniteBanded(lapack.HermitianUpper, OverlapCholeskyUpper);

	/*  From BLAS documentation
	 *  Lower hermitian banded storage
	 *
	 *  DO 20, J = 1, N
	 *      M = 1 - J
	 *      DO 10, I = J, MIN( N, J + K )
	 *          A( M + I, J ) = matrix( I, J )
	 *      ENDDO
	 *  ENDDO
  	 */
	OverlapHermitianLower.resize(BasisSize, SuperDiagonals+1);
	OverlapHermitianLower = 0.0;		
	for (int j=0; j<BasisSize; j++)
	{
		for (int i=j; i<std::min(BasisSize, j+SuperDiagonals+1); i++)
		{
			OverlapHermitianLower(j, i-j) = evaluator(i , j);
		}
	}


	//Set up "full" storage formats, used by InnerProduct
	int bandwidth = 2*SuperDiagonals + 1;
	int bw = SuperDiagonals + 1;
	OverlapFullCol.resize(BasisSize, bandwidth);
	OverlapFullRow.resize(bandwidth, BasisSize);
	OverlapFullCol = 0;
	OverlapFullRow = 0;
	for (int band=-bw+1; band<bw; band++)
	{
		int start=0;
		int end=BasisSize;
		if (band < 0) { start = -band; }
		else { end -= band; }
		for (int index=start; index<end; index++)
		{
			if (band<0)
			{
				int Jb = index + band;
				int Ib = - band;
				OverlapFullCol(index, band+bw-1) = real(OverlapHermitianLower(Jb, Ib)); 
				OverlapFullRow(band+bw-1, index) = real(OverlapHermitianLower(Jb, Ib));
			}
			else
			{
				int Jb = index;
				int Ib = band;
				OverlapFullCol(index, band+bw-1) = real(OverlapHermitianLower(Jb, Ib)); 
				OverlapFullRow(band+bw-1, index) = real(OverlapHermitianLower(Jb, Ib)); 
			}
		}
	}

}


/*
 * Vector
 */

void OverlapMatrix::MultiplyOverlapVector(const VectorType &source, VectorType &dest)
{
	BLAS<cplx> blas;
	blas.MultiplyMatrixVectorBandedHermitian(MatrixHermitianStorage::Upper, OverlapHermitianUpper, 1.0, source, 0.0, dest);
}

void OverlapMatrix::MultiplyOverlapVector(VectorType &vector)
{
	BLAS<cplx> blas;
	blas.MultiplyMatrixVectorBandedTriangular(MatrixHermitianStorage::Upper, MatrixTranspose::None, MatrixDiagonal::NonUnit, OverlapCholeskyUpper, vector);		
	blas.MultiplyMatrixVectorBandedTriangular(MatrixHermitianStorage::Upper, MatrixTranspose::Conjugate, MatrixDiagonal::NonUnit, OverlapCholeskyUpper, vector);
}

void OverlapMatrix::SolveOverlapVector(VectorType &vector)
{
	BLAS<cplx> blas;
	blas.SolveMatrixVectorBandedTriangular(MatrixHermitianStorage::Upper, MatrixTranspose::Conjugate, MatrixDiagonal::NonUnit, OverlapCholeskyUpper, vector);
	blas.SolveMatrixVectorBandedTriangular(MatrixHermitianStorage::Upper, MatrixTranspose::None, MatrixDiagonal::NonUnit, OverlapCholeskyUpper, vector);
}

void OverlapMatrix::MultiplySqrtOverlapVector(bool conj, VectorType &vector)
{
	BLAS<cplx> blas;
	MatrixTranspose transpose = conj ? MatrixTranspose::Conjugate : MatrixTranspose::None;
	blas.MultiplyMatrixVectorBandedTriangular(MatrixHermitianStorage::Upper, transpose, MatrixDiagonal::NonUnit, OverlapCholeskyUpper, vector);		
}

void OverlapMatrix::SolveSqrtOverlapVector(bool conj, VectorType &vector)
{
	BLAS<cplx> blas;
	MatrixTranspose transpose = conj ? MatrixTranspose::Conjugate : MatrixTranspose::None;
	blas.SolveMatrixVectorBandedTriangular(MatrixHermitianStorage::Upper, transpose, MatrixDiagonal::NonUnit, OverlapCholeskyUpper, vector);		
}


/*
 * Tensor
 */
void OverlapMatrix::MultiplyOverlapTensor(const TensorType &source, TensorType &dest)
{
	for (int i = 0; i < source.extent(0); i++)
	{
		for (int j = 0; j < source.extent(2); j++)
		{
			VectorType srcSlice = source(i, Range::all(), j);
			VectorType dstSlice = dest(i, Range::all(), j);
			MultiplyOverlapVector(srcSlice, dstSlice);
		}
	}
}

void OverlapMatrix::MultiplyOverlapTensor(TensorType &vector)
{
	for (int i = 0; i < vector.extent(0); i++)
	{
		for (int j = 0; j < vector.extent(2); j++)
		{
			Array<cplx, 1> slice = vector(i, Range::all(), j);
			MultiplyOverlapVector(slice);
		}
	}
}

void OverlapMatrix::SolveOverlapTensor(TensorType &vector)
{
	for (int i = 0; i < vector.extent(0); i++)
	{
		for (int j = 0; j < vector.extent(2); j++)
		{
			Array<cplx, 1> slice = vector(i, Range::all(), j);
			SolveOverlapVector(slice);
		}
	}
}

void OverlapMatrix::MultiplySqrtOverlapTensor(bool conj, TensorType &vector)
{
	for (int i = 0; i < vector.extent(0); i++)
	{
		for (int j = 0; j < vector.extent(2); j++)
		{
			Array<cplx, 1> slice = vector(i, Range::all(), j);
			MultiplySqrtOverlapVector(conj, slice);
		}
	}
}

void OverlapMatrix::SolveSqrtOverlapTensor(bool conj, TensorType &vector)
{
	for (int i = 0; i < vector.extent(0); i++)
	{
		for (int j = 0; j < vector.extent(2); j++)
		{
			Array<cplx, 1> slice = vector(i, Range::all(), j);
			SolveSqrtOverlapVector(conj, slice);
		}
	}
}




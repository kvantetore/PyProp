
#ifndef _BLITZLAPACK_H
#define _BLITZLAPACK_H

#include <blitz/array.h>
#include <complex>

#include "../common.h"
#include "fortran.h"

#define LAPACK_NAME(x) x ## _

namespace blitz
{
namespace linalg
{

extern "C"
{
	void LAPACK_NAME(zgeev)( char* JOBVL, char* JOBVR, int* N, cplx* A, int* LDA, cplx* W, cplx* VL, int* LDVL, cplx* VR, int* LDVR, cplx* WORK, int* LWORK, double* RWORK, int* INFO );
	void LAPACK_NAME(zheev)( char* JOBVZ, char* UPLO, int* N, cplx* A, int* LDA, double* W, cplx* WORK, int* LWORK, double* RWORK, int* INFO );
	void LAPACK_NAME(zgeqrf)( int* M, int* N, cplx* A, int* LDA, cplx* TAU, cplx* WORK, int* LWORK, int* INFO );
	void LAPACK_NAME(zungqr)( int* M, int* N, int *K, cplx* A, int* LDA, cplx* TAU, cplx* WORK, int* LWORK, int* INFO );
    void LAPACK_NAME(zunmqr)( char* SIDE, char* TRANS, int* M, int* N, int* K, cplx* A, int* LDA, cplx* TAU, cplx* C, int* LDC, cplx* WORK, int* LWORK, int* INFO );
	void LAPACK_NAME(zhbgv)( char* JOBZ, char* UPLO, int* N, int* KA, int* KB, cplx* AB, int* LDAB, cplx* BB, int* LDBB, double* W, cplx* Z, int* LDZ, cplx* WORK, double* RWORK, int* INFO);
	void LAPACK_NAME(zgbsv)(int* N, int* KL, int* KU, int* NRHS, cplx* AB, int* LDAB, int* IPIV, cplx* B, int* LDB, int* INFO);

	void LAPACK_NAME(zpbsvx)( char* FACT, char* UPLO, int* N, int* KD, int* NRHS, cplx* AB, int* LDAB, cplx* AFB, int* LDAFB, char* EQUED, double* S, cplx* B, int* LDB, cplx* X, int* LDX, double* RCOND, double* FERR, double* BERR, cplx* WORK, double* RWORK, int* INFO );

	void LAPACK_NAME(zpbtrf)(char* UPLO, int* N, int* KD, cplx* AB, int* LDAB, int* INFO );
	void LAPACK_NAME(zpbtrs)(char* UPLO, int* N, int* KD, int* NRHS, cplx* AB, int* LDAB, cplx* B, int* LDB, int* INFO );

	void LAPACK_NAME(zgetri)( int* N, cplx* A, int* LDA, int* IPIV, cplx* WORK, int* LWORK, int* INFO );
	void LAPACK_NAME(zgetrf)( int* N, int* M, cplx* A, int* LDA, int* IPIV, int* INFO );
	void LAPACK_NAME(zgetrs)( char* TRANS, int* N, int* NRHS, cplx* A, int* LDA, int* IPIV, cplx* B, int* LDB, int* INFO );

	void LAPACK_NAME(zgbtrs)( char* TRANS, int* N, int* KL, int* KU, int* NRHS, cplx* AB, int* LDAB, int* IPIV, cplx* B, int* LDB, int* INFO ); 
	void LAPACK_NAME(zgbtrf)( int* M, int* N, int* KL, int* KU, cplx* AB, int* LDAB, int* IPIV, int* INFO );
}


using namespace blitz;

template <class T>
class LAPACK
{

public:
	typedef Array<T, 1> VectorType;
	typedef Array<int, 1> VectorTypeInt;
	typedef Array<double, 1> DoubleVectorType;
	typedef Array<cplx, 1> ComplexDoubleVectorType;
	typedef Array<T, 2> MatrixType;
	typedef Array<double, 2> MatrixTypeDouble;

	LAPACK()
	{
		complexDoubleWork.resize(10);
		complexSingleWork.resize(10);
		doubleWork.resize(10);
		singleWork.resize(10);
	}

	enum MultiplicationSide
	{
		MultiplyLeft,
		MultiplyRight
	};

	enum MultiplicationTranspose
	{
		TransposeNone,
		TransposeConjugate
	};

	enum HermitianStorage
	{
		HermitianUpper,
		HermitianLower
	};

	enum MatrixFactoring
	{
		MatrixIsFactored,
		MatrixIsNotFactored,
		MatrixIsNotEquilibrated
	};

	enum MatrixEquilibrium
	{
		MatrixEquilibriumEnabled,
		MatrixEquilibriumDisabled
	};

	/*
	 * General Dense
	 */
	//LU Factorization
	int CalculateMatrixInverse(MatrixType &matrix, VectorTypeInt &pivot);
	int CalculateLUFactorization(MatrixType &matrix, VectorTypeInt &pivot);
	int SolveGeneralFactored(MultiplicationTranspose transpose, MatrixType &factoredMatrix, VectorTypeInt &pivot, MatrixType &rightHandSide);
	int SolveGeneralFactored(MultiplicationTranspose transpose, MatrixType &factoredMatrix, VectorTypeInt &pivot, VectorType &rightHandSide)
	{
		TinyVector<int, 2> shape(1, rightHandSide.extent(0));
		TinyVector<int, 2> stride(rightHandSide.extent(0)*rightHandSide.stride(0), rightHandSide.stride(0));
		MatrixType rhsMatrix(rightHandSide.data, shape, stride, neverDeleteData);
		return SolveGeneralFactored(transpose, factoredMatrix, pivot, rhsMatrix);
	}

	//Eigenvector Factorization
	int CalculateEigenvectorFactorization(bool calculateLeft, bool calculateRight, 
			MatrixType &matrix, VectorType &eigenvalues, MatrixType &leftEigenvectors, 
			MatrixType &rightEigenvectors);

	//QR Factorization
	int CalculateQRFactorization(MatrixType &matrix, VectorType &reflectors);
	int CompleteQRFactorization(MatrixType &matrix, VectorType &reflectors);
	int ApplyQRTransformation( MultiplicationSide side, MultiplicationTranspose transpose, MatrixType &qrMatrix, VectorType &reflectors, MatrixType &matrix );
	int ApplyQRTransformation( MultiplicationSide side, MultiplicationTranspose transpose, MatrixType &qrMatrix, VectorType &reflectors, VectorType &vector );

	
	/*
	 * Hermitian Dense
	 */
	int CalculateEigenvectorFactorizationHermitian(bool calculateEigenvectors, HermitianStorage storage, 
			MatrixType &matrix, DoubleVectorType &eigenvalues);


	/* 
	 * Hermitian Banded
	 */
	int CalculateEigenvectorFactorizationGeneralizedHermitianBanded(bool calculateEigenvetors, HermitianStorage storage, MatrixType &matrixA, 
			MatrixType &matrixB, DoubleVectorType &eigenvalues, MatrixType &eigenvectors);


	/*
	 * Positive Definite Banded
	 */
	int CalculateCholeskyFactorizationPositiveDefiniteBanded(HermitianStorage storage, MatrixType equationMatrix);
	int SolvePositiveDefiniteBandedFactored(HermitianStorage storage, MatrixType factoredMatrix, MatrixType rightHandSide);
	int SolvePositiveDefiniteBandedFactored(MultiplicationTranspose transpose, MatrixType &factoredMatrix, VectorTypeInt &pivot, VectorType &rightHandSide)
	{
		TinyVector<int, 2> shape(1, rightHandSide.extent(0));
		TinyVector<int, 2> stride(rightHandSide.extent(0)*rightHandSide.stride(0), rightHandSide.stride(0));

		MatrixType rhsMatrix(rightHandSide.data, shape, stride, neverDeleteData);
		return SolvePositiveDefiniteBandedFactored(transpose, factoredMatrix, pivot, rhsMatrix);
	}

	
	/* 
	 * General Banded 
	 */
	int SolveGeneralBandedSystemOfEquations(MatrixType &equationMatrix, VectorTypeInt &pivotVector, VectorType &rightHandVector, int &subDiagonals, int &superDiagonals);
	int CalculateLUFactorizationBanded(MatrixType equationMatrix, VectorTypeInt pivots);
	int SolveBandedFactored(MatrixType factoredMatrix, VectorTypeInt pivots, VectorType rightHandSide);


	/*
	 * Expert interface
	 */
	int SolvePositiveDefiniteBandedSystemOfEquations(MatrixFactoring fact, HermitianStorage storage, MatrixType equationMatrix, MatrixType factoredMatrix, char* EQUED, DoubleVectorType S, MatrixType rightHandSide, MatrixType solution, double* conditionEstimate, DoubleVectorType forwardErrorEstimate, DoubleVectorType backwardErrorEstimate);



private:
	Array<std::complex<double>, 1> complexDoubleWork;
	Array<std::complex<float>, 1> complexSingleWork;
	Array<double, 1> doubleWork;
	Array<float, 1> singleWork;

	void PreconditionSolveGeneralFactored(MultiplicationTranspose transpose, MatrixType &factoredMatrix, VectorTypeInt &pivot, MatrixType &rightHandSide);

	void PreconditionCalculateEigenvectorFactorization(bool calculateLeft, bool calculateRight, 
		MatrixType &matrix, VectorType &eigenvalues, MatrixType &leftEigenvectors, 
		MatrixType &rightEigenvectors);
	void PreconditionCalculateEigenvectorFactorizationHermitian(bool calculateEigenvectors, HermitianStorage storage, 
			MatrixType &matrix, DoubleVectorType &eigenvalues);
	
	void PreconditionCalculateQRFactorization(MatrixType &A, VectorType &reflectors);
	//void PreconditionCompleteQRFactorization(MatrixType R, VectorType reflectors, MatrixType Q);
	void PreconditionCompleteQRFactorization(MatrixType &matrix, VectorType &reflectors);
	void PreconditionApplyQRTransformation( MultiplicationSide side, MultiplicationTranspose transpose, MatrixType &qrMatrix, VectorType &reflectors, MatrixType &matrix );

	void PreconditionCalculateEigenvectorFactorizationGeneralizedHermitianBanded(bool calculateEigenvetors, HermitianStorage storage, MatrixType &matrixA, 
			MatrixType &matrixB, DoubleVectorType &eigenvalues, MatrixType &eigenvectors);

	void PreconditionSolveGeneralBandedSystemOfEquations(MatrixType &equationMatrix, VectorTypeInt &pivotVector, VectorType &rightHandVector, int &subDiagonals, int &superDiagonals);

	void PreconditionSolvePositiveDefiniteBandedSystemOfEquations(MatrixType equationMatrix, MatrixType factoredMatrix, DoubleVectorType S, MatrixType rightHandSide, MatrixType solution, DoubleVectorType ForwardErrorEstimate, DoubleVectorType BackwardErrorEstimate);
	
};


/* 
 * Check wheter all parameters are valid to the BLAS functions. 
 * This is type independent, and need only be written once
 */

template<class T>
void LAPACK<T>::PreconditionSolveGeneralFactored(MultiplicationTranspose transpose, MatrixType &factoredMatrix, VectorTypeInt &pivot, MatrixType &rightHandSide)
{
	BZPRECONDITION(factoredMatrix.extent(0) == factoredMatrix.extent(1));
	BZPRECONDITION(factoredMatrix.stride(1) == 1);
	BZPRECONDITION(pivot.extent(0) == factoredMatrix.extent(0));
	BZPRECONDITION(pivot.stride(0) == 1);
	BZPRECONDITION(rightHandSide.extent(1) == factoredMatrix.extent(0));
	BZPRECONDITION(rightHandSide.stride(1) == 1);
}

template<class T>
void LAPACK<T>::PreconditionCalculateEigenvectorFactorizationHermitian(bool calculateEigenvectors, HermitianStorage storage, 
			MatrixType &matrix, DoubleVectorType &eigenvalues)
{
	BZPRECONDITION(matrix.extent(0) == matrix.extent(1));
	BZPRECONDITION(matrix.extent(0) == eigenvalues.extent(0));
}
	

template<class T>
void LAPACK<T>::PreconditionCalculateEigenvectorFactorization(bool calculateLeft, bool calculateRight, 
		MatrixType &matrix, VectorType &eigenvalues, MatrixType &leftEigenvectors, 
		MatrixType &rightEigenvectors)
{
	BZPRECONDITION(matrix.extent(0) == matrix.extent(1));
	BZPRECONDITION(matrix.extent(0) == eigenvalues.extent(0));
	if (calculateLeft)
	{
		BZPRECONDITION(leftEigenvectors.extent(0) == leftEigenvectors.extent(1));
		BZPRECONDITION(leftEigenvectors.extent(0) == matrix.extent(0));
		BZPRECONDITION(leftEigenvectors.stride(1) == 1);
	}
	if (calculateRight)
	{
		BZPRECONDITION(rightEigenvectors.extent(0) == rightEigenvectors.extent(1));
		BZPRECONDITION(rightEigenvectors.extent(0) == matrix.extent(0));
		BZPRECONDITION(rightEigenvectors.stride(1) == 1);
	}

	BZPRECONDITION(matrix.stride(1) == 1);
	BZPRECONDITION(eigenvalues.stride(0) == 1);
}

template<class T>
void LAPACK<T>::PreconditionCalculateQRFactorization(MatrixType &matrix, VectorType &reflectors)
{
	BZPRECONDITION(std::min(matrix.extent(0), matrix.extent(0)) == reflectors.extent(0));
	BZPRECONDITION(matrix.stride(1) == 1);
}

template<class T>
void LAPACK<T>::PreconditionCompleteQRFactorization(MatrixType &matrix, VectorType &reflectors)
{
	BZPRECONDITION(std::min(matrix.extent(0), matrix.extent(0)) == reflectors.extent(0));
	BZPRECONDITION(matrix.stride(1) == 1);
}

template<class T>
void LAPACK<T>::PreconditionApplyQRTransformation( MultiplicationSide side, MultiplicationTranspose transpose, MatrixType &qrMatrix, VectorType &reflectors, MatrixType &matrix )
{
	//TODO: ADD preconditions
}


template<class T>
void LAPACK<T>::PreconditionCalculateEigenvectorFactorizationGeneralizedHermitianBanded(bool calculateEigenvetors, HermitianStorage storage, MatrixType &matrixA, 
			MatrixType &matrixB, DoubleVectorType &eigenvalues, MatrixType &eigenvectors)

{
	BZPRECONDITION(matrixA.extent(0) == matrixB.extent(0));
	BZPRECONDITION(matrixA.extent(1) == matrixB.extent(1));
	BZPRECONDITION(matrixA.extent(0) == eigenvalues.extent(0));
	BZPRECONDITION(matrixA.extent(0) == eigenvectors.extent(0));
}

template<class T>
void LAPACK<T>::PreconditionSolveGeneralBandedSystemOfEquations(MatrixType &equationMatrix, VectorTypeInt &pivotVector, VectorType &rightHandVector, 
		int &subDiagonals, int &superDiagonals)
{
	//TODO: Add preconditions
}

template<class T>
void LAPACK<T>::PreconditionSolvePositiveDefiniteBandedSystemOfEquations(MatrixType equationMatrix, MatrixType factoredMatrix, DoubleVectorType S, MatrixType rightHandSide, MatrixType solution, DoubleVectorType ForwardErrorEstimate, DoubleVectorType BackwardErrorEstimate)
{
	//Get data matrix sizes
//	int n = equationMatrix.extent(0);
//	int kd = equationMatrix.extent(1);
//	int ldab = kd;
//	int nrhs = rightHandSide.extent(0);
//	int ldafb = factoredMatrix.extent(1);
//	int ldb = rightHandSide.extent(1);
//	int ldx = solution.extent(1);

	BZPRECONDITION(factoredMatrix.extent(0) == equationMatrix.extent(0));
	BZPRECONDITION(equationMatrix.extent(1) == factoredMatrix.extent(1));
	
}

/*
 * Implementation for complex<double>
 */


template<>
inline int LAPACK<cplx>::CalculateLUFactorization(MatrixType &matrix, blitz::Array<int, 1> &pivot)
{
	//TODO: add preconditioner

	int N = matrix.extent(0);
	int M = matrix.extent(1);
	int LDA = M;

	int info = 0;

	//Call LAPACK routine to compute LU factorization
	LAPACK_NAME(zgetrf)(&N, &M, matrix.data(), &LDA, pivot.data(), &info);

	if (info != 0)
	{
		cout << "WARNING: zgetrf could not compute LU factorization, info = " << info << endl;
	}

	return info;
}

template<>
inline int LAPACK<cplx>::SolveGeneralFactored(MultiplicationTranspose transpose, MatrixType &factoredMatrix, VectorTypeInt &pivot, MatrixType &rightHandSide)
{
	PreconditionSolveGeneralFactored(transpose, factoredMatrix, pivot, rightHandSide);

	char transChar = transpose == TransposeNone ? 'N' : 'C';
	int N = factoredMatrix.extent(1);
	int lda = factoredMatrix.stride(0);
	int ldb = rightHandSide.stride(0);
	int nrhs = rightHandSide.extent(0);
	int info;

	LAPACK_NAME(zgetrs)( &transChar, &N, &nrhs, factoredMatrix.data(), &lda, pivot.data(), rightHandSide.data(), &ldb, &info );
	if (info != 0)
	{
		cout << "WARNING: zgetrs could not back substitute LU factorization, info = " << info << endl;
	}	
	return info;
}


template<>
inline int LAPACK<cplx>::CalculateMatrixInverse(MatrixType &matrix, blitz::Array<int, 1> &pivot)
{
	//TODO: add preconditioner
	
	int N = matrix.extent(0);
	int LDA = matrix.extent(1);

	int info = 0;
	int workLength = N;

	//Resize work arrays if needed
	if (complexDoubleWork.extent(0) < N)
	{
		complexDoubleWork.resize(N);
	}

	//Call LAPACK routine to invert matrix
	LAPACK_NAME(zgetri)(&N, matrix.data(), &LDA, pivot.data(), complexDoubleWork.data(), &workLength, &info);

	if (info != 0)
	{
		cout << "WARNING: zgetri could not invert matrix, info = " << info << endl; 
	}

	return info;
}


template<>
inline int LAPACK<cplx>::CalculateEigenvectorFactorizationHermitian(bool calculateEigenvectors, HermitianStorage storage, 
			MatrixType &matrix, DoubleVectorType &eigenvalues)
{
	PreconditionCalculateEigenvectorFactorizationHermitian(calculateEigenvectors, storage, matrix, eigenvalues);

	char jobz = calculateEigenvectors ? 'V' : 'N';
	char uplo = storage == HermitianUpper ? 'U' : 'L';

	int N = matrix.extent(0);
	int LDA = matrix.stride(0);

	int info;
	int workLength = -1;

	//Make sure realwork is big enough
	if (doubleWork.size() < std::max(1, 3*N-2))
	{
		doubleWork.resize(std::max(1, 3*N-2));
	}

	//First call: calculate optimal work array size
	LAPACK_NAME(zheev)(&jobz, &uplo, &N, matrix.data(), &LDA, eigenvalues.data(), 
		complexDoubleWork.data(), &workLength, doubleWork.data(), &info );

	//Resize work arrays
	int optimalWorkLength = (int) real(complexDoubleWork(0));
	if (complexDoubleWork.size() < optimalWorkLength)
	{
		complexDoubleWork.resize(optimalWorkLength);
	}

	//Second call: calculate eigenvalues/vectors
	workLength = complexDoubleWork.extent(0);//optimalWorkLength;
	LAPACK_NAME(zheev)(&jobz, &uplo, &N, matrix.data(), &LDA, eigenvalues.data(), 
		complexDoubleWork.data(), &workLength, doubleWork.data(), &info );


	if (info != 0)
	{
		cout << "WARNING: zgeev could not find requested eigenvalues: " << info << endl;
	}
	return info;
}


template<>
inline int LAPACK<cplx>::CalculateEigenvectorFactorization(bool calculateLeft, bool calculateRight, 
		MatrixType &matrix, VectorType &eigenvalues, MatrixType &leftEigenvectors, 
		MatrixType &rightEigenvectors)
{
	PreconditionCalculateEigenvectorFactorization(calculateLeft, calculateRight, matrix, 
		eigenvalues, leftEigenvectors, rightEigenvectors);


	char leftChar = calculateLeft ? 'V' : 'N';
	char rightChar = calculateRight ? 'V' : 'N';

	int N = matrix.extent(0);
	int LDA = matrix.stride(0);
	int LDL = leftEigenvectors.stride(0);
	int LDR = rightEigenvectors.stride(0);

	int info;
	int workLength = -1;

	//First call: calculate optimal work array size
	LAPACK_NAME(zgeev)(&leftChar, &rightChar, &N, matrix.data(), &LDA, eigenvalues.data(), 
		leftEigenvectors.data(), &LDL, rightEigenvectors.data(), &LDR, complexDoubleWork.data(), &workLength, 
		doubleWork.data(), &info);

	//Resize work arrays
	int optimalWorkLength = (int) real(complexDoubleWork(0));
	if (doubleWork.size() < optimalWorkLength)
	{
		doubleWork.resize(optimalWorkLength);
	}
	if (complexDoubleWork.size() < optimalWorkLength)
	{
		complexDoubleWork.resize(optimalWorkLength);
	}

	//Second call: calculate eigenvalues/vectors
	workLength = complexDoubleWork.extent(0);//optimalWorkLength;
	LAPACK_NAME(zgeev)(
		&leftChar, 
		&rightChar, 
		&N, 
		matrix.data(), 
		&LDA, 
		eigenvalues.data(), 
		leftEigenvectors.data(), 
		&LDL, 
		rightEigenvectors.data(), 
		&LDR, 
		complexDoubleWork.data(), 
		&workLength, 
		doubleWork.data(), 
		&info);

	if (info != 0)
	{
		cout << "WARNING: zgeev could not find requested eigenvalues: " << info << endl;
	}
	return info;
}


template<>
inline int LAPACK<cplx>::CalculateQRFactorization(MatrixType &matrix, VectorType &reflectors)
{
	PreconditionCalculateQRFactorization(matrix, reflectors);

	//matrix is in row-major (FORTRAN) format
	int m = matrix.extent(1);
	int n = matrix.extent(0);
	int lda = matrix.stride(0);
	int info;

	//First run, find optimal work array length
	int workLength = -1;
	LAPACK_NAME(zgeqrf)( &m, &n, matrix.data(), &lda, reflectors.data(), complexDoubleWork.data(), &workLength, &info );

	int optimalWorkLength = (int) real(complexDoubleWork(0));
	optimalWorkLength = std::max(optimalWorkLength, n) ;
	if (complexDoubleWork.size() < optimalWorkLength)
	{
		complexDoubleWork.resize(optimalWorkLength);
	}

	//Second run, calculate!
	workLength = optimalWorkLength;
	LAPACK_NAME(zgeqrf)( &m, &n, matrix.data(), &lda, reflectors.data(), complexDoubleWork.data(), &workLength, &info );

	if (info != 0)
	{
		cout << "WARNING: Calculation of QR factorization failed: " << info << endl;
	}

	return info;
}

template<>
inline int LAPACK<cplx>::CompleteQRFactorization(MatrixType &matrix, VectorType &reflectors)
{
	PreconditionCompleteQRFactorization(matrix, reflectors);

	//matrix is in row-major (FORTRAN) format
	int m = matrix.extent(1);
	int n = matrix.extent(0);
	int k = std::min(m,n);
	int lda = matrix.stride(0);
	int info;

	//First run, find optimal work array length
	int workLength = -1;
	LAPACK_NAME(zungqr)( &m, &n, &k, matrix.data(), &lda, reflectors.data(), complexDoubleWork.data(), &workLength, &info );

	int optimalWorkLength = (int) real(complexDoubleWork(0));
	optimalWorkLength = std::max(optimalWorkLength, n) ;
	if (complexDoubleWork.size() < optimalWorkLength)
	{
		complexDoubleWork.resize(optimalWorkLength);
	}

	//Second run, calculate!
	workLength = optimalWorkLength;
	LAPACK_NAME(zungqr)( &m, &n, &k, matrix.data(), &lda, reflectors.data(), complexDoubleWork.data(), &workLength, &info );

	if (info != 0)
	{
		cout << "WARNING: Calculation of QR factorization failed: " << info << endl;
	}

	return info;
}


template<>
inline int LAPACK<cplx>::ApplyQRTransformation( MultiplicationSide side, MultiplicationTranspose transpose, MatrixType &qrMatrix, VectorType &reflectors, MatrixType &matrix )
{
	PreconditionApplyQRTransformation(side, transpose, qrMatrix, reflectors, matrix);

	//matrix is in row-major (FORTRAN) format
	int m = matrix.extent(1);
	int n = matrix.extent(0);
	int k = std::min(qrMatrix.extent(0), qrMatrix.extent(1));
	int lda = qrMatrix.stride(0);
	int ldc = matrix.stride(0);
	int info;

	char sideChar = side == MultiplyLeft ? 'L' : 'R';
	char transChar = transpose == TransposeNone ? 'N' : 'C';

	//First run, calculate worklength
	int workLength = -1;

    LAPACK_NAME(zunmqr)( &sideChar, &transChar, &m, &n, &k, qrMatrix.data(), &lda, reflectors.data(), matrix.data(), &ldc, complexDoubleWork.data(), &workLength, &info );

	int optimalWorkLength = (int) real(complexDoubleWork(0));
	optimalWorkLength = std::max(optimalWorkLength, n) ;
	if (complexDoubleWork.size() < optimalWorkLength)
	{
		complexDoubleWork.resize(optimalWorkLength);
	}

	//Second run, calculate!
	workLength = optimalWorkLength;
    LAPACK_NAME(zunmqr)( &sideChar, &transChar, &m, &n, &k, qrMatrix.data(), &lda, reflectors.data(), matrix.data(), &ldc, complexDoubleWork.data(), &workLength, &info );

	return info;
}

template<>
inline int LAPACK<cplx>::ApplyQRTransformation( MultiplicationSide side, MultiplicationTranspose transpose, MatrixType &qrMatrix, VectorType &reflectors, VectorType &vector )
{
	blitz::TinyVector<int, 2> shape(1, vector.extent(0));
	blitz::TinyVector<int, 2> stride(vector.extent(0)*vector.stride(0), vector.stride(0));

	MatrixType matrix(vector.data(), shape, stride, blitz::neverDeleteData);
	return ApplyQRTransformation(side, transpose, qrMatrix, reflectors, matrix);
}


/*
 * Solve generalized eigenvalue problem for hermitian banded matrices A and B:  Ax = lBx
 */
template<>
inline int LAPACK<cplx>::CalculateEigenvectorFactorizationGeneralizedHermitianBanded(bool calculateEigenvectors,
                                                                                     HermitianStorage storage,
                                                                                     MatrixType &matrixA, 
                                                                                     MatrixType &matrixB, 
                                                                                     DoubleVectorType &eigenvalues, 
                                                                                     MatrixType &eigenvectors)
{

	// Some sanity checks on input
	PreconditionCalculateEigenvectorFactorizationGeneralizedHermitianBanded(calculateEigenvectors, storage, matrixA, matrixB, eigenvalues, eigenvectors);

	// Should we calculate eigenvectors
	char jobz = calculateEigenvectors ? 'V' : 'N';

	// How is matrices stored?
	char uplo = storage == HermitianUpper ? 'U' : 'L';

	// Dimensions of input matrices
	int n = matrixA.extent(0);
	int ldab = matrixA.extent(1);
	int ldbb = matrixB.extent(1);
	int ka = ldab - 1;
	int kb = ldbb - 1;

	int ldz = n;

	// Resize work arrays if needed
	if (complexDoubleWork.extent(0) < n) { complexDoubleWork.resize(n); }
	if (doubleWork.extent(0) < 3 * n) { doubleWork.resize(3 * n); }

	int info;

	
	//Call LAPACK to perfrom calculation
	LAPACK_NAME(zhbgv)(
		&jobz,                       // Eigenvalues only or eigenvalues+eigenvectors
		&uplo,                       // What is stored: upper or lower triangle of matrix
		&n,                          // Order of input matrices >= 0
		&ka,                         // Number of super- / sub-diagonals in matrix A
		&kb,                         // Number of super- / sub-diagonals in matrix B
		matrixA.data(),
		&ldab,                       // Leading dimension of A, >= KA+1
		matrixB.data(),
		&ldbb,                       // Leading dimension of B, >= KB+1
		eigenvalues.data(),          // (out) Eigenvalues in ascending order
		eigenvectors.data(),         // (out) Eigenvectors (Z)
		&ldz,                        // (in) Leading dimension of eigenvetor array
		complexDoubleWork.data(),
		doubleWork.data(),
		&info);


	if (info != 0)
	{
		if (info < 0)
		{
			cout << "WARNING: (zhbgv) argument i had illegal value, i = " << info << endl;
		}
		else if (info <= n)
		{
			cout << "WARNING: (zhbgv) algorithm failed to converge! INFO =  " << info << endl;
		}
		else if (info > n)
		{
			cout << "WARNING: (zhbgv) matrix B not positive definite! INFO = " << info << endl;
		}
	}	 

	return info;
}

/*
 * Solve banded linear system of complex equations.
 */
template<>
inline int LAPACK<cplx>::SolveGeneralBandedSystemOfEquations(MatrixType &equationMatrix, VectorTypeInt &pivotVector, VectorType &rightHandVector, int &subDiagonals, int &superDiagonals)
{
	// Dimensions of input arrays
	int n = equationMatrix.extent(0);
	int kl = subDiagonals;
	int ku = superDiagonals;
	int nhrs = 1;
	int ldab = equationMatrix.extent(1);
	int ldb = rightHandVector.extent(0);

	int info;
	
	// Call LAPACK routine, solve system of equations
	LAPACK_NAME(zgbsv) ( &n, &kl, &ku, &nhrs, equationMatrix.data(), &ldab, pivotVector.data(), rightHandVector.data(), &ldb, &info);

	if (info != 0)
	{
		cout << "WARNING: Could not solve system of equations " << info << endl;
	}

	return info;
}

/*
 * Solve positive definite banded system of equations using LAPACK expert interface
 */
template<>
inline int LAPACK<cplx>::SolvePositiveDefiniteBandedSystemOfEquations(MatrixFactoring fact, 
                                                                      HermitianStorage storage, 
                                                                      MatrixType equationMatrix, 
                                                                      MatrixType factoredMatrix, 
                                                                      char* EQUED, 
																	  DoubleVectorType S, 
                                                                      MatrixType rightHandSide, 
                                                                      MatrixType solution, 
                                                                      double* conditionEstimate, 
                                                                      DoubleVectorType ForwardErrorEstimate, 
                                                                      DoubleVectorType BackwardErrorEstimate)
{

	PreconditionSolvePositiveDefiniteBandedSystemOfEquations(equationMatrix, factoredMatrix, S, rightHandSide, solution, 
		ForwardErrorEstimate, BackwardErrorEstimate);

	char FACT = (fact == MatrixIsFactored) ? 'F' : ( (fact == MatrixIsNotEquilibrated) ? 'E' : 'N');
	char uplo = storage == HermitianUpper ? 'U' : 'L';

	//Get data matrix sizes
	int n = equationMatrix.extent(0);
	int ldab = equationMatrix.extent(1);
	int kd = ldab - 1;
	int nrhs = rightHandSide.extent(0);
	int ldafb = factoredMatrix.extent(1);
	int ldb = rightHandSide.extent(1);
	int ldx = solution.extent(1);
	int info;

	//Resize work arrays if needed
	if (complexDoubleWork.size() < 2 * n)
	{
		complexDoubleWork.resize(2 * n);
	}

	if (doubleWork.size() < n)
	{
		doubleWork.resize(n);
	}

	LAPACK_NAME(zpbsvx)(
		&FACT, 
		&uplo, 
		&n, 
		&kd, 
		&nrhs, 
		equationMatrix.data(), 
		&ldab, 
		factoredMatrix.data(), 
		&ldafb, 
		EQUED, 
		S.data(), 
		rightHandSide.data(), 
		&ldb, 
		solution.data(), 
		&ldx, 
		conditionEstimate, 
		ForwardErrorEstimate.data(), 
		BackwardErrorEstimate.data(), 
		complexDoubleWork.data(),
		doubleWork.data(), 
		&info 
		);


		if (info != 0)
		{
			if (info < 0)
			{
				cout << "WARNING: (zpbsvx) argument i had illegal value, i = " << info << endl;
			}
			else if (info <= n)
			{
				cout << "WARNING: (zpbsvx) algorithm failed to converge! INFO =  " << info << endl;
			}
			else if (info > n)
			{
				cout << "WARNING: (zpbsvx) matrix singular to working precision! INFO = " << info << endl;
			}
		}

	return info;
}

template<>
inline int LAPACK<cplx>::CalculateCholeskyFactorizationPositiveDefiniteBanded(HermitianStorage storage, MatrixType equationMatrix)
{
	char uplo = storage == HermitianUpper ? 'U' : 'L';

	//Get data matrix sizes
	int n = equationMatrix.extent(0);
	int ldab = equationMatrix.extent(1);
	int kd = ldab - 1;
	int info;

	LAPACK_NAME(zpbtrf)(&uplo, &n, &kd, equationMatrix.data(), &ldab, &info );

	if (info != 0)
	{
		cout << "WARNING: zpbrtf error " << info << endl;
	}

	return info;
}

template<>
inline int LAPACK<cplx>::SolvePositiveDefiniteBandedFactored(HermitianStorage storage, MatrixType factoredMatrix, MatrixType rightHandSide)
{
	char uplo = storage == HermitianUpper ? 'U' : 'L';

	//Get data matrix sizes
	int n = factoredMatrix.extent(0);
	int ldab = factoredMatrix.extent(1);
	int kd = ldab - 1;
	int nrhs = rightHandSide.extent(0);
	int ldb = rightHandSide.extent(1);
	int info;

	LAPACK_NAME(zpbtrs)(&uplo, &n, &kd, &nrhs, factoredMatrix.data(), &ldab, rightHandSide.data(), &ldb, &info );
	if (info != 0)
	{
		cout << "WARNING: zpbrts error " << info << endl;
	}

	return info;
}

template<>
inline int LAPACK<cplx>::CalculateLUFactorizationBanded(MatrixType equationMatrix, VectorTypeInt pivots)
{
	/* 
	 * Assumptions (IMPORTANT)
	 *   - equationMatrix is square
	 *   - number of superdiagonal bands == number of subdiagonal bands
	 */
	
	//Get data matrix sizes
	int m = equationMatrix.extent(0);
	int n = m;
	int ldab = equationMatrix.extent(1);
	int kl = (ldab - 1) / 3;
	int ku = kl;

	int info;

	LAPACK_NAME(zgbtrf)(&m, &n, &kl, &ku, equationMatrix.data(), &ldab, pivots.data(), &info);

	if (info != 0)
	{
		cout << "WARNING: zgbrtf error " << info << endl;
	}

	return info;
}

template<>
inline int LAPACK<cplx>::SolveBandedFactored(MatrixType factoredMatrix, VectorTypeInt pivots, VectorType rightHandSide)
{
	//Get data matrix sizes
	int m = factoredMatrix.extent(0);
	int n = m;
	int ldab = factoredMatrix.extent(1);
	int kl = (ldab - 1) / 3;
	int ku = kl;
	int nrhs = 1;
	int ldb = rightHandSide.size();

	char trans = 'N';
	int info;

	LAPACK_NAME(zgbtrs)(&trans, &n, &kl, &ku, &nrhs, factoredMatrix.data(), &ldab, pivots.data(), rightHandSide.data(), &ldb, &info);

	if (info != 0)
	{
		cout << "WARNING: zgbrts error " << info << endl;
	}

	return info;
}


}} //Namespaces

#endif


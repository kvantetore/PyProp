#ifndef PIRAM_BLITZLAPACK_H
#define PIRAM_BLITZLAPACK_H

#include <blitz/array.h>
#include <complex>

#include <core/common.h>

#include "blitzblas.h"


#define LAPACK_NAME(x) x ## _

namespace blitz
{
namespace linalg
{

extern "C"
{
	void LAPACK_NAME(zgeev)( char* JOBVL, char* JOBVR, int* N, cplx* A, int* LDA, cplx* W, cplx* VL, int* LDVL, cplx* VR, int* LDVR, cplx* WORK, int* LWORK, double* RWORK, int* INFO );
	void LAPACK_NAME(zgeqrf)( int* M, int* N, cplx* A, int* LDA, cplx* TAU, cplx* WORK, int* LWORK, int* INFO );
	void LAPACK_NAME(zungqr)( int* M, int* N, int *K, cplx* A, int* LDA, cplx* TAU, cplx* WORK, int* LWORK, int* INFO );
    void LAPACK_NAME(zunmqr)( char* SIDE, char* TRANS, int* M, int* N, int* K, cplx* A, int* LDA, cplx* TAU, cplx* C, int* LDC, cplx* WORK, int* LWORK, int* INFO );
}


using namespace blitz;

template <class T>
class LAPACK
{

public:
	typedef Array<T, 1> VectorType;
	typedef Array<T, 2> MatrixType;

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

	int CalculateEigenvectorFactorization(bool calculateLeft, bool calculateRight, 
			MatrixType &matrix, VectorType &eigenvalues, MatrixType &leftEigenvectors, 
			MatrixType &rightEigenvectors);

	int CalculateQRFactorization(MatrixType &matrix, VectorType &reflectors);
	//void CompleteQRFactorization(MatrixType R, VectorType reflectors, MatrixType Q);
	int CompleteQRFactorization(MatrixType &matrix, VectorType &reflectors);
	int ApplyQRTransformation( MultiplicationSide side, MultiplicationTranspose transpose, MatrixType &qrMatrix, VectorType &reflectors, MatrixType &matrix );
	int ApplyQRTransformation( MultiplicationSide side, MultiplicationTranspose transpose, MatrixType &qrMatrix, VectorType &reflectors, VectorType &vector );

private:
	Array<std::complex<double>, 1> complexDoubleWork;
	Array<std::complex<float>, 1> complexSingleWork;
	Array<double, 1> doubleWork;
	Array<float, 1> singleWork;
	BLAS<T> blas;

	void PreconditionCalculateEigenvectorFactorization(bool calculateLeft, bool calculateRight, 
		MatrixType &matrix, VectorType &eigenvalues, MatrixType &leftEigenvectors, 
		MatrixType &rightEigenvectors);
	
	void PreconditionCalculateQRFactorization(MatrixType &A, VectorType &reflectors);
	//void PreconditionCompleteQRFactorization(MatrixType R, VectorType reflectors, MatrixType Q);
	void PreconditionCompleteQRFactorization(MatrixType &matrix, VectorType &reflectors);
	void PreconditionApplyQRTransformation( MultiplicationSide side, MultiplicationTranspose transpose, MatrixType &qrMatrix, VectorType &reflectors, MatrixType &matrix );
};


/* 
 * Check wheter all parameters are valid to the BLAS functions. 
 * This is type independent, and need only be written once
 */

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

/*
 * Implementation for complex<double>
 */

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

}} //Namespaces

#endif


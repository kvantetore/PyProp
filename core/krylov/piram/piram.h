#ifndef PIRAM_H
#define PIRAM_H

#include <core/common.h>


template <class T, class OP>
class pIRAM
{
public:
	typedef blitz::Array<T, 1> VectorType;
	typedef blitz::Array<int, 1> IntVectorType;
	typedef blitz::Array<T, 2> MatrixType;

private:
	OP matrixOperator;

	//Parameters
	int MaxRestartCount;
	int MaxOrthogonalizationCount;
	int EigenvalueCount;
	int MatrixSize
	int BasisSize;
	int Tolerance; 
	//...

	//Iteration variables
	/*
	 * Matrices in pIRAM are always stored in col-major transposed
	 * format (transposed C style). For a matrix M, this means
	 * M(col, row) will give you the expected element of M. Furthermore,
	 * it means that M(col, row+1) will be the consecutive element in
	 * memory, making it more cache efficient, and easier to interoperate
	 * with fortran libraries like LAPACK and BLAS
	 */
	//Large matrices/vectors (~ size of MatrixSize)
	MatrixType ArnoldiVectors;
	VectorType Residual;
	VectorType TempVector;
	//Small matrices/vectors (~ size of BasisSize)
	MatrixType HessenbergMatrix;
	MatrixType HessenbergEigenvectors;
	MatrixType HessenbergTriangular;
	VectorType Overlap;
	VectorType Overlap2;
	MatrixType Eigenvalues;
	IntVectorType EigenvalueOrdering;
	//Scalars
	int CurrentArnoldiStep;
	bool IsConverged;
	//...
	
	//Statistics
	int OperatorCount;
	int OrthogonalizationCount;
	int RestartCount;
	

public:
	pIRAM();

	void Setup();
	void SetupResidual();
	void MultiplyOperator(VectorType &in, VectorType &out);
	void PerformArnoldiStep();
	

};

template <class T, class OP>
void pIRAM<T, OP>::Setup()
{
	ArnoldiVectors.resize(MatrixSize, BasisSize);
	HessenbergMatrix.resize(BasisSize, BasisSize);
	Residual.resize(MatrixSize);
	TempVector.resize(MatrixSize);
	HessenbergEigenvectors.resize(BasisSize, BasisSize);
	HessenbergTriangular.resize(BasisSize, BasisSize);
	Overlap.resize(BasisSize);
	Overlap2.resize(BasisSize);
	Eigenvalues.resize(BasisSize);
	EigenvalueOrdering.resize(BasisSize);

	ArnoldiVectors = 0;
	HessenbergMatrix = 0;
	CurrentArnoldiStep = 0;
	IsConverged = false;

	//Perform initial arnoldi step
	SetupResidual();
	VectorType v0(ArnoldiVectors(0, blitz::Range::all()));
	CopyVector(Residual, v0);
	MultiplyOperator(v0, Residual);
	T alpha = InnerProduct(v0, Residual);
	AddVector(v0, -alpha, Residual);
	self.HessenbergMatrix(0,0) = alpha;
}

/*
 * CopyVector
 * InnerProduct
 * AddVector
 */

template <class T, class OP>
void pIRAM<T, OP>::PerformArnoldiStep()
{
	//Define some shortcuts
	int j = CurrentArnoldiStep;
	VectorType currentArnoldiVector(ArnoldiVectors(j+1, blitz::Range::all()));
	MatrixType currentArnoldiMatrix(ArnoldiVectors(blitz::Range(0, j+1), blitz::Range::all()));
	VectorType currentOverlap(Overlap(0, j+1));

	//Update the Arnoldi Factorization with the previous Residual
	T beta = VectorNorm(Residual);
	ScaleVector(Residual, norm);
	CopyVector(Residual, currentArnoldiVector);
	HessenbergMatrix(j, j+1) = beta;

	//Expand the krylov supspace with 1 dimension
	MultiplyOperator(currentArnoldiVector, Residual);
	T origNorm = VectorNorm(Residual);

	//Perform Gram Schmidt to remove components of Residual 
	//which are linear combinations of the j first Arnoldi vectors
	//- Calculate the projection of Residual into the range of currentArnoldiMatrix (B)
	//  i.e. TempVector =  B* B Residual
	MultiplyMatrixVector(MatrixTranspose::Conjugate, currentArnoldiMatrix, Residual, currentOverlap);
	MultiplyMatrixVector(MatrixTranspose::None, currentArnoldiMatrix, currentOverlap, TempVector);
	//- Remove the projection from the residual
	AddVector(TempVector, -1.0, Residual);
	T residualNorm = VectorNorm(Residual);

	//If necessary, perform any reorthogonalization steps to ensure that all 
	//arnoldi vectors are orthogonal
	PerformOrthogonalization(origNorm, residualNorm)
	
	//Update the Hessenberg Matrix
	//TODO: Use CopyVector
	HessenbergMatrix(j+1, blitz::Range(0, j+1)) = currentOverlap;

	CurrentArnoldiStep++;
}

template <class T, class OP>
void pIRAM<T, OP>::PerformOrthogonalization(T originalNorm, T residualNorm)
{
	//Define some shortcuts
	int j = CurrentArnoldiStep;
	VectorType currentArnoldiVector(ArnoldiVectors(j+1, blitz::Range::all()));
	MatrixType currentArnoldiMatrix(ArnoldiVectors(blitz::Range(0, j+1), blitz::Range::all()));
	VectorType currentOverlap(Overlap(0, j+1));
	VectorType currentOverlap2(Overlap2(0, j+1));

	/*
	 * 0.717 ~= sin( pi/4 )
	 */
	int curOrthoCount = 0;
	while (residualNorm < 0.717 * originalNorm)
	{
		//Update statistics
		OrthogonalizationCount++;
		curOrthoCount++;

		if (curOrthoCount > MaxOrthogonalizationCount)
		{
			cout << "WARNING: Maximum orthogonalization steps for one Arnoldi iteration has been reached" << endl;
			cout << "         Spurious eigenvalues may occur in the solution" << endl;
			break;
		}

		//Perform one step gram schmidt step to keep vectors orthogonal
		MultiplyMatrixVector(MatrixTranspose::Conjugate, currentArnoldiMatrix, Residual, currentOverlap2);
		MultiplyMatrixVector(MatrixTranspose::None, currentArnoldiMatrix, currentOverlap2, TempVector);
		//- Remove the projection from the residual
		AddVector(TempVector, -1.0, Residual);
		AddVector(currentOverlap2, 1.0, currentOverlap);
	
		//Prepare for next iteration
		originalNorm = residualNorm;
		residualNorm = VectorNorm(Residual);		
	}
}


template <class T, class OP>
void UpdateEigenvalues()
{
	//Calculate eigenvalues with the eigenvector factorization
	CalculateEigenvectorFactorization(HessenbergMatrix, HessenbergEigenvectors, Eigenvalues);
	SortVector(Eigenvalues, EigenvaluesOrdering);

	/*
	 * Check if solution is converged. See ARPACK Users Guide for information
	 * about the convergence cirterium
	 */
	double tol = Tolerance;
	double eps = GetMachinePrecision<double>()
	double normHessenberg = ...;
	double normResidual = VectorNorm(Residual);

	IsConverged = true;
	for (int i=0; i<EigenvalueCount; i++)
	{
		int sortedIndex = EigenvaluesOrdering(i);
		double ritzError = HessenbergEigenvectors(sortedIndex, BasisSize);
		double ritzValue = Eigenvalues(sortedIndex);
		double convergenceCrit = fabs(ritzError) * normResidual - std::max(eps * normHessenberg, tol * abs(ritzValue));

		IsConverged = IsConverged && convergenceCrit<0;
	}
}

template <class T, class OP>
void PerformRestartStep()
{
	//Update statistics
	RestartCount++;

	VectorType q(BasisSize);
	q = 0.0;
	q(BasisSize-1) = 1.0;

	//Use the p=BasisSize - EigenvalueCount eigenvalues of the Hessenberg matrix
	//as shifts in a truncated shifted QR algorithm, in order to compress the 
	//interesting eigenvectors into the BasisSize first arnoldi vectors
	for (int i=EigenvalueCount; i<BasisSize; i++)
	{
		//Get shift
		T shift = Eigenvalues(EigenvalueOrdering(i));

		//Shift the Hessenberg matrix
		for (int j=0; j<BasisSize; j++)
		{
			HessenbergMatrix(j,j) -= shift;
		}

		//TODO: Implement this efficiantly for Hessenberg matrices 
		//(by using bulge chase givens rotations or similar)
		
		//Calculate QR factorization of the shifted matrix
		//reuse the HessenbergEigenvectors as basisMatrix for the QR 
		MatrixType basisMatrix(HessenbergEigenvectors);
		CalculateQR(HessenbergMatrix, basisMatrix, HessenbergTriangular);

		//Apply the unitary transformation on HessenbergMatrix, ArnoldiVectors and q
		MatrixMatrixMultiplication(HessenbergTriangular, basisMatrix, HessenbergMatrix)
		//TODO: Make sure this works in-place
		MatrixMatrixMultiplication(ArnoldiVectors, basisMatrix, ArnoldiVectors)
		MatrixVectorMultiplication(MatrixTranspose::Conjugate, basisMatrix, q);
		q = conj(q);
	}

	//Restart the arnoldi iteration on step k
	int k = EigenvalueCount-1;
	self.CurrentArnoldiStep = k;

	//Calculate residual for our restarted arnoldi iteration
	VectorType arnoldiVector(ArnoldiVectors(k+1, blitz::Range::all()));
	T sigma = q(k);
	T beta = HessenbergMatrix(k, k+1);
	ScaleVector(Residual, sigma);
	AddVector(arnoldiVector, beta, Residual);

	//Clear the last p cols in HessenbergMatrix
	HessenbergMatrix( blitz::Range(k+1, BasisSize), blitz::Range::all() ) = 0;
}

#endif


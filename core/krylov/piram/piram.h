#ifndef PIRAM_H
#define PIRAM_H

#include <limits>
#include <core/common.h>
#include <core/mpi/mpitraits.h>
#include <random/uniform.h>
#include "blitzblas.h"
#include "blitzlapack.h"

#include <iostream>
#include <fstream>
#include <strstream>


namespace piram
{

using namespace blitz::linalg;

template <class T, class OP>
class pIRAM
{
public:
	typedef blitz::Array<T, 1> VectorType;
	typedef blitz::Array<int, 1> IntVectorType;
	typedef blitz::Array<T, 2> MatrixType;
	typedef blitz::linalg::LAPACK<T> LAPACK;
	typedef blitz::linalg::BLAS<T> BLAS;

	//Parameters
	int MaxRestartCount;
	int MaxOrthogonalizationCount;
	int EigenvalueCount;
	int MatrixSize;
	int BasisSize;
	double Tolerance;
	MPI_Comm CommBase;
	//...

	//Constructor
	pIRAM() 
	{
		//TODO: Implement default values
	}

private:
	OP matrixOperator;
	LAPACK lapack;
	BLAS blas;


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
	VectorType Overlap3;
	VectorType Reflectors;
	VectorType Eigenvalues;
	IntVectorType EigenvalueOrdering;
	VectorType ErrorEstimates;
	//Scalars
	int CurrentArnoldiStep;
	bool IsConverged;
	//...
	
	//Statistics
	int OperatorCount;
	int OrthogonalizationCount;
	int RestartCount;

	//Private methods
	void PerformOrthogonalization(T originalNorm, T residualNorm);
	void UpdateEigenvalues();
	void PerformRestartStep();
	void SortVector(VectorType &vector, IntVectorType &ordering);
	void RecreateArnoldiFactorization();
	T CalculateGlobalNorm(VectorType &vector)
	{
		return sqrt(CalculateGlobalInnerProduct(vector, vector));
	}
	
	T CalculateGlobalInnerProduct(VectorType &x, VectorType &y)
	{
		//T localNorm = blas.VectorNorm(vector);
		T localValue = blas.InnerProduct(x, y);
		T globalValue;

		MPITraits<T> traits;
		int status = MPI_Allreduce(&localValue, &globalValue, 1*traits.Length(), traits.Type(), MPI_SUM, MPI_COMM_WORLD);

		return globalValue;
	}

	void GlobalAddVector(VectorType &in, VectorType &out)
	{
		MPITraits<T> traits;
		int status = MPI_Allreduce(in.data(), out.data(), in.extent(0)*traits.Length(), traits.Type(), MPI_SUM, MPI_COMM_WORLD);
	}

public:
	void Solve();
	void Setup();
	void SetupResidual();
	void MultiplyOperator(VectorType &in, VectorType &out);
	void PerformArnoldiStep();
};

template <class T, class OP>
void pIRAM<T, OP>::Setup()
{
	ArnoldiVectors.resize(BasisSize, MatrixSize);
	Residual.resize(MatrixSize);
	TempVector.resize(MatrixSize);
	HessenbergMatrix.resize(BasisSize, BasisSize);
	HessenbergEigenvectors.resize(BasisSize, BasisSize);
	HessenbergTriangular.resize(BasisSize, BasisSize);
	Overlap.resize(BasisSize);
	Overlap2.resize(BasisSize);
	Overlap3.resize(BasisSize);
	Eigenvalues.resize(BasisSize);
	EigenvalueOrdering.resize(BasisSize);
	Reflectors.resize(BasisSize);
	ErrorEstimates.resize(EigenvalueCount);

	//Reset statistcs
	OperatorCount = 0;
	RestartCount = 0;
	OrthogonalizationCount = 0;

	//Reset arnoldi factorization
	ArnoldiVectors = 0;
	HessenbergMatrix = 0;
	CurrentArnoldiStep = 0;
	IsConverged = false;

	//Perform initial arnoldi step
	SetupResidual();
	VectorType v0(ArnoldiVectors(0, blitz::Range::all()));
	blas.CopyVector(Residual, v0);
	MultiplyOperator(v0, Residual);
	T alpha = CalculateGlobalInnerProduct(v0, Residual);
	blas.AddVector(v0, -alpha, Residual);
	HessenbergMatrix(0,0) = alpha;
}

template <class T, class OP>
void pIRAM<T, OP>::PerformArnoldiStep()
{
	//Define some shortcuts
	int j = CurrentArnoldiStep;
	VectorType currentArnoldiVector(ArnoldiVectors(j+1, blitz::Range::all()));
	MatrixType currentArnoldiMatrix(ArnoldiVectors(blitz::Range(0, j+1), blitz::Range::all()));
	VectorType currentOverlap(Overlap(blitz::Range(0, j+1)));
	VectorType currentOverlap2(Overlap2(blitz::Range(0, j+1)));

	//Update the Arnoldi Factorization with the previous Residual
	T beta = CalculateGlobalNorm(Residual);
	blas.ScaleVector(Residual, 1.0/beta);
	blas.CopyVector(Residual, currentArnoldiVector);
	HessenbergMatrix(j, j+1) = beta;

	//Expand the krylov supspace with 1 dimension
	MultiplyOperator(currentArnoldiVector, Residual);
	T origNorm = CalculateGlobalNorm(Residual);

	//Perform Gram Schmidt to remove components of Residual 
	//which are linear combinations of the j first Arnoldi vectors
	//- Calculate the projection of Residual into the range of currentArnoldiMatrix (B)
	//  i.e. TempVector =  B* B Residual
	blas.MultiplyMatrixVector(MatrixTranspose::Conjugate, currentArnoldiMatrix, 1.0, Residual, 0.0, currentOverlap2);
	GlobalAddVector(currentOverlap2, currentOverlap);
	
	blas.MultiplyMatrixVector(currentArnoldiMatrix, currentOverlap, TempVector);
	//- Remove the projection from the residual
	blas.AddVector(TempVector, -1.0, Residual);
	T residualNorm = CalculateGlobalNorm(Residual);

	//If necessary, perform any reorthogonalization steps to ensure that all 
	//arnoldi vectors are orthogonal
	PerformOrthogonalization(origNorm, residualNorm);
	
	//Update the Hessenberg Matrix
	//TODO: Use CopyVector
	HessenbergMatrix(j+1, blitz::Range(0, j+1)) = currentOverlap;

	/*
	int procId, procCount;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);
	for (int j=0; j<procCount; j++)
	{
		if (procId == j)
		{
			cout << "ProcId = " << procId << ", Hessenberg = " << real(HessenbergMatrix) << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	throw std::runtime_error("hei");
	*/

	CurrentArnoldiStep++;
}

template <class T, class OP>
void pIRAM<T, OP>::PerformOrthogonalization(T originalNorm, T residualNorm)
{
	//Define some shortcuts
	int j = CurrentArnoldiStep;
	VectorType currentArnoldiVector(ArnoldiVectors(j+1, blitz::Range::all()));
	MatrixType currentArnoldiMatrix(ArnoldiVectors(blitz::Range(0, j+1), blitz::Range::all()));
	VectorType currentOverlap(Overlap(blitz::Range(0, j+1)));
	VectorType currentOverlap2(Overlap2(blitz::Range(0, j+1)));
	VectorType currentOverlap3(Overlap3(blitz::Range(0, j+1)));

	/*
	 * 0.717 ~= sin( pi/4 )
	 */
	int curOrthoCount = 0;
	while (std::abs(residualNorm) < 0.717 * std::abs(originalNorm))
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
		blas.MultiplyMatrixVector(MatrixTranspose::Conjugate, currentArnoldiMatrix, 1.0, Residual, 0.0, currentOverlap2);
		GlobalAddVector(currentOverlap2, currentOverlap3);
		blas.MultiplyMatrixVector(currentArnoldiMatrix, currentOverlap3, TempVector);
		//- Remove the projection from the residual
		blas.AddVector(TempVector, -1.0, Residual);
		blas.AddVector(currentOverlap3, 1.0, currentOverlap);
	
		//Prepare for next iteration
		originalNorm = residualNorm;
		residualNorm = CalculateGlobalNorm(Residual);		
	}
}


template <class T, class OP>
void pIRAM<T, OP>::SortVector(VectorType &vector, IntVectorType &ordering)
{
	BZPRECONDITION(vector.extent(0) == ordering.extent(0));
	ordering = blitz::tensor::i;

	//TODO: For f*** sake, don't implement our own sort algo.  

	int n = vector.extent(0);
	for (int i=0; i<n; i++)
	{
		for (int j=i+1; j<n; j++)
		{
			if (std::abs(vector(i)) > std::abs(vector(j)))
			{
				T temp = vector(i);
				vector(i) = vector(j);
				vector(j) = temp;

				int itemp = ordering(i);
				ordering(i) = ordering(j);
				ordering(j) = itemp;
			}
		}
	}
}


template <class T, class OP>
void pIRAM<T,OP>::UpdateEigenvalues()
{
	MatrixType emtpyMatrix(1,1);


	//Calculate eigenvalues with the eigenvector factorization
	//CalculateEigenvectorFactorization destroys the input matrix, so we make a copy
	HessenbergTriangular = HessenbergMatrix; 
	lapack.CalculateEigenvectorFactorization(false, true, HessenbergTriangular, Eigenvalues, emtpyMatrix, HessenbergEigenvectors);
	SortVector(Eigenvalues, EigenvalueOrdering); 

	/*
	 * Check if solution is converged. See ARPACK Users Guide for information
	 * about the convergence cirterium
	 */
	double tol = Tolerance;
	double eps = std::numeric_limits<double>::epsilon();
	double normHessenberg = 0; //TODO: calculate some norm here
	double normResidual = std::abs(CalculateGlobalNorm(Residual));

	IsConverged = true;
	for (int i=0; i<EigenvalueCount; i++)
	{
		int sortedIndex = EigenvalueOrdering(i);
		double ritzError = std::abs(HessenbergEigenvectors(sortedIndex, BasisSize-1));
		double ritzValue = std::abs(Eigenvalues(i));
		double convergenceCrit = ritzError * normResidual - std::max(eps * normHessenberg , tol * ritzValue);
		ErrorEstimates(i) = convergenceCrit;
		IsConverged = IsConverged && convergenceCrit<0;
	}
}

template <class T, class OP>
void pIRAM<T,OP>::PerformRestartStep()
{
	//Update statistics
	RestartCount++;

	//We need q for the restarted residual
	VectorType q(BasisSize);
	q = 0.0;
	q(BasisSize-1) = 1.0;

	//Use the p=BasisSize - EigenvalueCount eigenvalues of the Hessenberg matrix
	//as shifts in a truncated shifted QR algorithm, in order to compress the 
	//interesting eigenvectors into the BasisSize first arnoldi vectors
	for (int i=EigenvalueCount; i<BasisSize; i++)
	{
		//Get shift
		T shift = Eigenvalues(i);

		//Shift the Hessenberg matrix
		for (int j=0; j<BasisSize; j++)
		{
			HessenbergMatrix(j,j) -= shift;
		}

		//TODO: Implement this efficiantly for Hessenberg matrices 
		//(by using bulge chase givens rotations or similar)
		
		//Calculate QR factorization of the shifted matrix
		lapack.CalculateQRFactorization(HessenbergMatrix, Reflectors);
		HessenbergTriangular = HessenbergMatrix;

		//the upper triangular part of HessenbergMatrix is now R, while
		//the lower triangular part are the vectors v forming Q by HH reflections

		//Apply the unitary transformation on HessenbergMatrix, ArnoldiVectors and q
		lapack.ApplyQRTransformation(LAPACK::MultiplyRight, LAPACK::TransposeNone, HessenbergTriangular, Reflectors, ArnoldiVectors);
		lapack.ApplyQRTransformation(LAPACK::MultiplyLeft, LAPACK::TransposeConjugate, HessenbergTriangular, Reflectors, q);
		q = conj(q);

		//Clear the subdiagonals on HessenbergMatrix
		for (int j=0; j<BasisSize-1; j++)
		{
			HessenbergMatrix(j, blitz::Range(j+1, BasisSize-1)) = 0;
		}
		//Left-multiply with Q
		lapack.ApplyQRTransformation(LAPACK::MultiplyRight, LAPACK::TransposeNone, HessenbergTriangular, Reflectors, HessenbergMatrix);
		for (int j=0; j<BasisSize; j++)
		{
			HessenbergMatrix(j, j) += shift;
		}
	
	}

//		cout << "Hessenberg = " << real(HessenbergMatrix) << endl;
	//Restart the arnoldi iteration on step k
	int k = EigenvalueCount-1;
	CurrentArnoldiStep = k;

	//Calculate residual for our restarted arnoldi iteration
	VectorType arnoldiVector(ArnoldiVectors(k+1, blitz::Range::all()));
	T sigma = q(k);
	T beta = HessenbergMatrix(k, k+1);
	blas.ScaleVector(Residual, sigma);
	blas.AddVector(arnoldiVector, beta, Residual);
	//Clear the last p cols in HessenbergMatrix to make it hessenberg
	HessenbergMatrix( blitz::Range(k+1, BasisSize-1), blitz::Range::all() ) = 0;
}


template<class T, class OP>
void pIRAM<T, OP>::MultiplyOperator(VectorType &in, VectorType &out)
{
	//Update statistics
	OperatorCount++;

	//Perform matrix-vector multiplication
	matrixOperator(in, out);
}


template<class T, class OP>
void pIRAM<T, OP>::SetupResidual()
{
	using namespace ranlib;
	using namespace std;

	UniformClosed<double> rand;

	int procId, procCount;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);

	int localStart = procId * MatrixSize;
	cout << "proc " << procId <<": LocalStart = " << localStart << endl;

	ifstream file;
	file.open("init.dat", ios::in);

	int i=0;
	int localIndex = 0;
	while (! file.eof() && file.is_open())
	{
		
		char line[256];
		file.getline(line, 256);

		if (i>=localStart && i<localStart+MatrixSize)
		{
			std::istrstream lineStream(line);
			lineStream >> Residual(localIndex);;
			localIndex++;
		}
		
		i++;
	}

	/*
	cout << "i = " << i << ", localIndex = " << localIndex << endl;

	for (int j=0; j<procCount; j++)
	{
		if (procId == j)
		{
			cout << "ProcId = " << procId << ", Residual = " << real(Residual) << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}


	if (procCount == 1)
	{
		VectorType v1 = Residual(blitz::Range(0, MatrixSize/2-1));
		VectorType v2 = Residual(blitz::Range(MatrixSize/2, MatrixSize));
		cout << "firstNorm = " << CalculateGlobalNorm(v1) << ", secondNorm = " << CalculateGlobalNorm(v2) << endl;
	}
	*/

	T norm = CalculateGlobalNorm(Residual);
	blas.ScaleVector(Residual, 1.0/norm);
	//throw std::runtime_error("hei");
	
	file.close();
}

template<class T, class OP>
void pIRAM<T, OP>::RecreateArnoldiFactorization()
{
	//Make sure we have a BasisSize-step arnoldi iteration
	while (CurrentArnoldiStep < BasisSize-1)
	{
		PerformArnoldiStep();
	}
}


template<class T, class OP>
void pIRAM<T, OP>::Solve()
{
	RecreateArnoldiFactorization();

	while (true)
	{
		UpdateEigenvalues();

		if (IsConverged)
		{
			int procId, procCount;
			MPI_Comm_rank(MPI_COMM_WORLD, &procId);
			MPI_Comm_size(MPI_COMM_WORLD, &procCount);

			if (procId == 0)
			{
				cout << "Converged!" << endl;
				cout << "Eigenvalues  = " << Eigenvalues(blitz::Range(0, EigenvalueCount-1)) << endl;
				cout << "OpCount      = " << OperatorCount << endl;
				cout << "RestartCount = " << RestartCount << endl;
				cout << "ReOrthoCOunt = " << OrthogonalizationCount << endl;
			}
			break;
		}

		PerformRestartStep();
		RecreateArnoldiFactorization();
	}
}




/*
 * BLAS:
 * CopyVector
 * VectorNorm
 * MatrixNorm //Check with ARPACK
 * InnerProduct
 * AddVector
 * ScaleVector
 * MultiplyMatrixVector //Check wheter these two can be performed in-place
 * MultiplyMatrixMatrix //
 *
 *
 * LAPACK:
 * CalculateEigenvectorFactorization
 * CalculateQR
 *
 * OTHER:
 * GetMachinePrecision
 * SortVector
 */

} //Namespace


#endif


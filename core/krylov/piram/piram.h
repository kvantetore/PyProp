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
	//Helper types
	typedef double NormType;
	typedef blitz::Array<T, 1> VectorType;
	typedef blitz::Array<NormType, 1> NormVectorType;
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
	NormType Tolerance;
	bool UseRandomStart;

	//Options
	MPI_Comm CommBase;
	bool DisableMPI;

	//Operator class, must implement operator()(VectorType &in, VectorType &out)
	OP matrixOperator;

	//Constructor
	pIRAM() 
	{
		Tolerance = std::numeric_limits<NormType>::epsilon();
		CommBase = MPI_COMM_WORLD;
		DisableMPI = false;
		UseRandomStart;
	}

private:
	LAPACK lapack;
	BLAS blas;

	int ProcId;
	int ProcCount;

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
	NormVectorType ErrorEstimates;
	//Scalars
	int CurrentArnoldiStep;
	bool IsConverged;
	//...
	
	//Statistics
	int OperatorCount;
	int OrthogonalizationCount;
	int RestartCount;

	//Private methods
	void MultiplyOperator(VectorType &in, VectorType &out);
	void PerformArnoldiStep();
	void PerformOrthogonalization(T originalNorm, T residualNorm);
	void UpdateEigenvalues();
	void PerformRestartStep();
	void SortVector(VectorType &vector, IntVectorType &ordering);
	void RecreateArnoldiFactorization();
	//Post-processing:
	void CalculateEigenvectors()
	{
		/* 
		 * If A V = V H is a an arnoldi factorization of A, 
		 * and H = X L inv(X) is an eigenvector factorization of H,
		 * A V X = V X L => Y = V X => A Y = Y L is an approximation 
		 * to a truncated eigenvector factorization of A
		 */

		/* 
		 * In order to save memory, we calculate Y = V X row-wise, and put the
		 * result back in ArnoldiVectors
		 */
		for (int i=0; i<MatrixSize; i++)
		{
			VectorType arnoldiRow = ArnoldiVectors(blitz::Range(0, BasisSize-1), i);
			//Use overlap as temporary vector and calculate the i-th row of Y
			blas.MultiplyMatrixVector(MatrixTranspose::Transpose, HessenbergEigenvectors, 1.0, arnoldiRow, 0.0, Overlap);
			arnoldiRow = Overlap;
		}
	}

	T CalculateGlobalNorm(VectorType &vector)
	{
		return sqrt(CalculateGlobalInnerProduct(vector, vector));
	}
	
	T CalculateGlobalInnerProduct(VectorType &x, VectorType &y)
	{
		T localValue = blas.InnerProduct(x, y);
		if (DisableMPI)
		{
			return localValue;
		}

		T globalValue;
		MPITraits<T> traits;
		MPI_Allreduce(&localValue, &globalValue, 1*traits.Length(), traits.Type(), MPI_SUM, CommBase);

		return globalValue;
	}

	void GlobalAddVector(VectorType &in, VectorType &out)
	{
		if (DisableMPI)
		{
			out = in;
		}

		MPITraits<T> traits;
		MPI_Allreduce(in.data(), out.data(), in.extent(0)*traits.Length(), traits.Type(), MPI_SUM, CommBase);
	}

public:
	void Setup();
	void Solve();
	void Postprocess()
	{
		CalculateEigenvectors();
	}

	void SetupResidual();

	//Accessor functions:
	int GetConvergedEigenvaluesCount()
	{
		return count(ErrorEstimates <= 1.0);
	}
	
	VectorType GetEigenvalues()
	{
		return Eigenvalues(blitz::Range(0, EigenvalueCount-1));
	}

	VectorType GetEigenvector(int eigenvectorIndex)
	{
		int sortedIndex = EigenvalueOrdering(eigenvectorIndex);
		return ArnoldiVectors(sortedIndex, blitz::Range(0, MatrixSize-1));
	}

	int GetOrthogonalizationCount()
	{
		return OrthogonalizationCount;
	}

	int GetRestartCount()
	{
		return RestartCount;
	}

	int GetOperatorCount()
	{
		return OperatorCount;
	}
		
};

template <class T, class OP>
void pIRAM<T, OP>::Setup()
{
	//Setup MPI
	if (!DisableMPI)
	{	
		MPI_Comm_rank(CommBase, &ProcId);
		MPI_Comm_size(CommBase, &ProcCount);
	}
	else
	{
		ProcId = 0;
		ProcCount = 1;
	}

	//Allocate workspace-memory
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
	NormType tol = Tolerance;
	NormType normResidual = std::abs(CalculateGlobalNorm(Residual));
	NormType eps = std::numeric_limits<NormType>::epsilon();
	NormType eps23 = std::pow(eps, 2.0/3.0);

	IsConverged = true;
	for (int i=0; i<EigenvalueCount; i++)
	{
		int sortedIndex = EigenvalueOrdering(i);
		NormType ritzError = std::abs(HessenbergEigenvectors(sortedIndex, BasisSize-1));
		NormType ritzValue = std::abs(Eigenvalues(i));
		NormType errorBounds = ritzError * normResidual;
		NormType accuracy = tol * std::max(eps23, tol * ritzValue);
		ErrorEstimates(i) = errorBounds / accuracy;
		IsConverged = IsConverged && (ErrorEstimates(i) <= 1);
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
	MPI_Comm_rank(CommBase, &procId);
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

	T norm = CalculateGlobalNorm(Residual);
	blas.ScaleVector(Residual, 1.0/norm);
	
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
			break;
		}

		if (RestartCount >= MaxRestartCount)
		{
			break;
		}

		PerformRestartStep();
		RecreateArnoldiFactorization();
	}


	if (ProcId == 0)
	{
		if (IsConverged)
		{
			cout << "All eigenvalues are converged!" << endl;
		}
		else
		{
			cout << "WARNING: Not all eigenvalues are converged" << endl;
			blitz::Array<bool, 1> conv (blitz::where(ErrorEstimates <= 1.0, true, false));
			cout << "    Converged Eigenvalues   = " << conv << endl;
			cout << "    Error Estimates (< 1.0) = " << ErrorEstimates << endl;
		}
		cout << "    Eigenvalues  = " << GetEigenvalues() << endl;
		cout << "    OpCount      = " << OperatorCount << endl;
		cout << "    RestartCount = " << RestartCount << endl;
		cout << "    ReOrthoCOunt = " << OrthogonalizationCount << endl;
	}
}

} //Namespace


#endif


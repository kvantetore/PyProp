#ifndef GMRES_H
#define GMRES_H

#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include <core/common.h>
#include <core/mpi/mpitraits.h>
#include <core/utility/blitzlapack.h>
#include <core/utility/timer.h>

#include "../piram/piram/blitzblas.h"
#include "../piram/piram/functors.h"


namespace gmres
{
using namespace blitz::linalg;

using namespace piram;

/*
 * Main GMRES-class
 */
template <class T>
class GMRES
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

	typedef std::map< std::string, Timer > TimerMap;

	//Parameters
	int MaxOrthogonalizationCount;
	int MatrixSize;
	int BasisSize;
	double Tolerance;

	bool PerformDoubleOrthogonalization;
	
	//Options
	MPI_Comm CommBase;
	bool DisableMPI;

	//Operator class, must implement operator()(VectorType &in, VectorType &out)
	typename OperatorFunctor<T>::Ptr MatrixOperator;
	typename IntegrationFunctor<T, NormType>::Ptr Integration;

	//Constructor
	GMRES() 
	{	
		//Default Params
		MaxOrthogonalizationCount = 3;
		CommBase = MPI_COMM_WORLD;
		DisableMPI = false;

		Tolerance = sqrt(std::numeric_limits<NormType>::epsilon());
		PerformDoubleOrthogonalization = true;
	}

private:
	LAPACK lapack;
	BLAS blas;

	int ProcId;
	int ProcCount;

	//Iteration variables
	/*
	 * Matrices in GMRES are always stored in col-major transposed
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
	MatrixType HessenbergSolve;
	VectorType HessenbergSolution;
	VectorType HessenbergSolution2;
	VectorType ErrorEstimateList;

	VectorType Overlap;
	VectorType Overlap2;
	VectorType Overlap3;
	//Scalars
	int CurrentArnoldiStep;
	bool IsConverged;
	bool HappyBreakdown;
	//...
	
	//Statistics
	int OperatorCount;
	int OrthogonalizationCount;
	int InnerProductCount;
	TimerMap Timers;

	//Error estimates
	T ErrorEstimate;

	//Private methods
	void MultiplyOperator(VectorType &in, VectorType &out);
	void PerformInitialArnoldiStep();
	void PerformArnoldiStep();
	void PerformOrthogonalization(T originalNorm, T residualNorm);
	void RecreateArnoldiFactorization();

	//Independent methods
	
	T CalculateGlobalNorm(VectorType &vector)
	{
		return Integration->Norm(vector);
	}
	
	T CalculateGlobalInnerProduct(VectorType &x, VectorType &y)
	{
		//Update Statistics
		InnerProductCount++;

		Timers["InnerProduct"].Start();
		T globalValue = Integration->InnerProduct(x, y);
		Timers["InnerProduct"].Stop();

		return globalValue;
	}


public:
	double EstimateMemoryUsage(int matrixSize, int basisSize);
	void Setup();
	void SolveVector(VectorType rhs, VectorType solution);
	void PrintStatistics();

	int GetOrthogonalizationCount()
	{
		return OrthogonalizationCount;
	}

	int GetInnerProductCount()
	{
		return InnerProductCount;
	}

	int GetOperatorCount()
	{
		return OperatorCount;
	}

	void ResetStatistics()
	{
		//Reset statistcs
		OperatorCount = 0;
		InnerProductCount = 0;
		OrthogonalizationCount = 0;

		Timers = TimerMap();
	}

	double GetResidualNorm()
	{
		return real(CalculateGlobalNorm(Residual));
	}

	MatrixType GetHessenbergMatrix()
	{
		return HessenbergMatrix;
	}

	double GetErrorEstimate() 
	{
		return abs(ErrorEstimate);
	}

	VectorType GetErrorEstimateList()
	{
		return ErrorEstimateList;
	}

		
};

template <class T>
double GMRES<T>::EstimateMemoryUsage(int matrixSize, int basisSize)
{
	double largeSize = matrixSize * (basisSize + 3.0);
	double smallSize = basisSize * (3*basisSize + 7.0);
	return (largeSize + smallSize) * sizeof(T) / (1024.0*1024.0);
}


template <class T>
void GMRES<T>::Setup()
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

	double minTol = std::pow(std::numeric_limits<NormType>::epsilon(), 2.0/3.0);
	if (Tolerance < minTol)
	{
		Tolerance = minTol;
	}

	//Use blas for InnerProducts unless another
	//integration functor has ben specified
	if (Integration == 0)
	{
		Integration = typename IntegrationFunctor<T, NormType>::Ptr(new BLASIntegrationFunctor<T, NormType>(DisableMPI, CommBase));
	}

	//Allocate workspace-memory
	ArnoldiVectors.resize(BasisSize, MatrixSize);
	Residual.resize(MatrixSize);
	TempVector.resize(MatrixSize);
	HessenbergMatrix.resize(BasisSize, BasisSize+1);
	HessenbergSolve.resize(BasisSize, BasisSize+1);
	HessenbergSolution.resize(BasisSize+1);
	HessenbergSolution2.resize(BasisSize+1);
	ErrorEstimateList.resize(BasisSize);
	Overlap.resize(BasisSize);
	Overlap2.resize(BasisSize);
	Overlap3.resize(BasisSize);

	ResetStatistics();

	//Reset arnoldi factorization
	ArnoldiVectors = 0;
	HessenbergMatrix = 0;
	CurrentArnoldiStep = 0;
}


template <class T>
void GMRES<T>::PerformInitialArnoldiStep()
{
	HappyBreakdown = false;

	ArnoldiVectors = 0; //This shouldn't be necessary
	HessenbergMatrix = 0;
	CurrentArnoldiStep = 0;

	//Normalize Residual
	//T norm = CalculateGlobalNorm(Residual);
	//cout << "Initial norm = " << norm << endl;
	//blas.ScaleVector(Residual, 1.0/norm);

	//First Arnoldi Vector is the init residual
	VectorType v0(ArnoldiVectors(0, blitz::Range::all()));
	blas.CopyVector(Residual, v0);

	//Apply operator
	MultiplyOperator(v0, Residual);

	//Update Hessenberg Matrix
	T alpha = CalculateGlobalInnerProduct(v0, Residual);
	blas.AddVector(v0, -alpha, Residual);
	HessenbergMatrix(0,0) = alpha;

	T beta = CalculateGlobalNorm(Residual);
	if (std::abs(beta) < Tolerance)
	{
		//cout << "Happy breakdown!, res = " << beta << ", alpha = " << alpha << endl;
		HappyBreakdown = true;
	}
	HessenbergMatrix(0, 1) = beta;

}

template <class T>
void GMRES<T>::PerformArnoldiStep()
{
	Timers["Arnoldi Step"].Start();

	//Define some shortcuts
	int j = CurrentArnoldiStep;
	VectorType currentArnoldiVector(ArnoldiVectors(j+1, blitz::Range::all()));
	MatrixType currentArnoldiMatrix(ArnoldiVectors(blitz::Range(0, j+1), blitz::Range::all()));
	VectorType currentOverlap(Overlap(blitz::Range(0, j+1)));
	VectorType currentOverlap2(Overlap2(blitz::Range(0, j+1)));

	//The previous normalized residual is the current start vector
	blas.CopyVector(Residual, currentArnoldiVector);
	T norm = HessenbergMatrix(j, j+1);
	blas.ScaleVector(currentArnoldiVector, 1.0/norm);

	//Expand the krylov supspace with 1 dimension
	Timers["Arnoldi Step"].Stop();
	MultiplyOperator(currentArnoldiVector, Residual);
	Timers["Arnoldi Step"].Start();
	T origNorm = CalculateGlobalNorm(Residual);

	//Perform Gram Schmidt to remove components of Residual 
	//which are linear combinations of the j first Arnoldi vectors
	//- Calculate the projection of Residual into the range of currentArnoldiMatrix (B)
	//  i.e. TempVector =  B* B Residual
	// Put the result in currentOverlap, use currentOverlap2 as temp buffer
	Timers["Arnoldi Step (Orthogonalization)"].Start();
	Integration->InnerProduct(currentArnoldiMatrix, Residual, currentOverlap, currentOverlap2);
	blas.MultiplyMatrixVector(currentArnoldiMatrix, currentOverlap, TempVector);
	//- Remove the projection from the residual
	blas.AddVector(TempVector, -1.0, Residual);
	T residualNorm = CalculateGlobalNorm(Residual);
	Timers["Arnoldi Step (Orthogonalization)"].Stop();

	//Perform additional orthogonalization (double orthogonalization)
	if (PerformDoubleOrthogonalization)
	{
		Timers["Arnoldi Step"].Stop();
		PerformOrthogonalization(origNorm, residualNorm);
		Timers["Arnoldi Step"].Start();
	}

	//Update the Hessenberg Matrix
	HessenbergMatrix(j+1, blitz::Range(0, j+1)) = currentOverlap;

	T beta = CalculateGlobalNorm(Residual);
	if (std::abs(beta) < Tolerance)
	{
		//cout << "Happy breakdown (2)!, res = " << beta << endl;
		HappyBreakdown = true;
	}
	HessenbergMatrix(j+1, j+2) = beta;

	CurrentArnoldiStep++;

	Timers["Arnoldi Step"].Stop();
}

template <class T>
void GMRES<T>::PerformOrthogonalization(T originalNorm, T residualNorm)
{
	Timers["Orthogonalization"].Start();

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
			cout << "WARNING: Maximum orthogonalization steps for one Arnoldi iteration has been " << endl
				 << "         reached (" << MaxOrthogonalizationCount << "), "
				 << "Spurious eigenvalues may occur in the solution" << endl;
			break;
		}

		//Perform one step gram schmidt step to keep vectors orthogonal
		// Put the result in currentOverlap2, use currentOverlap2 as temp buffer
		Integration->InnerProduct(currentArnoldiMatrix, Residual, currentOverlap3, currentOverlap2);
		blas.MultiplyMatrixVector(currentArnoldiMatrix, currentOverlap3, TempVector);
		//- Remove the projection from the residual
		blas.AddVector(TempVector, -1.0, Residual);
		blas.AddVector(currentOverlap3, 1.0, currentOverlap);
	
		//Prepare for next iteration
		originalNorm = residualNorm;
		residualNorm = CalculateGlobalNorm(Residual);		

		//cout << "ArnoldiStep = " << CurrentArnoldiStep << ", Residual Norm = " << residualNorm << endl;

	}

	Timers["Orthogonalization"].Stop();
}


template<class T>
void GMRES<T>::MultiplyOperator(VectorType &in, VectorType &out)
{
	Timers["Operator"].Start();

	//Update statistics
	OperatorCount++;

	//Perform matrix-vector multiplication
	(*MatrixOperator)(in, out);

	Timers["Operator"].Stop();
}


template<class T>
void GMRES<T>::RecreateArnoldiFactorization()
{
	//Make sure we have a BasisSize-step arnoldi iteration
	while (CurrentArnoldiStep < BasisSize-1 && !HappyBreakdown)
	{
		PerformArnoldiStep();
	}
}


template<class T>
void GMRES<T>::SolveVector(VectorType rightHandSide, VectorType solution)
{
	Timers["Total"].Start();
	Residual = rightHandSide;
	double inputNorm = std::abs(CalculateGlobalNorm(Residual));
	blas.ScaleVector(Residual, 1.0/inputNorm);

	/*
	 * A x = b
	 *
	 * x_n = Q_n y
	 * 
	 * min : || A Q_n y_n  - b ||
	 *
	 * min : || Q_n+1 H^_n y - b ||
	 *
	 * min : || H^_n y - ||b|| e_1 ||
	 *
	 */

	PerformInitialArnoldiStep();
	if (HappyBreakdown)
	{
		cplx c = HessenbergMatrix(0,0);
		HessenbergSolution(0) = inputNorm/c;
	}
	ErrorEstimate = HessenbergMatrix(CurrentArnoldiStep, CurrentArnoldiStep+1) ;
	ErrorEstimateList(CurrentArnoldiStep) = ErrorEstimate;

	while (CurrentArnoldiStep < BasisSize-1 && !HappyBreakdown && (std::abs(ErrorEstimate)>Tolerance))
	{
		PerformArnoldiStep();

		//Calculate minimal residual of current matrix
		Timers["Minimize Residual"].Start();
		int krylovSize = CurrentArnoldiStep + 1;
		HessenbergSolve = HessenbergMatrix;
		MatrixType H = HessenbergSolve( blitz::Range(0, krylovSize-1), blitz::Range(0, krylovSize) );
		VectorType y = HessenbergSolution( blitz::Range(0, krylovSize) );
		y = 0;
		y(0) = inputNorm;
		lapack.SolveLeastSquareGeneral(lapack.TransposeNone, H, y);
		VectorType x = HessenbergSolution( blitz::Range(0, krylovSize-1) );
		Timers["Minimize Residual"].Stop();
		
		//Calculate error
		Timers["Error Estimation"].Start();
		HessenbergSolve = HessenbergMatrix;
		VectorType r = HessenbergSolution2( blitz::Range(0, krylovSize) );
		blas.MultiplyMatrixVector(H, x, r);
		r(0) -= inputNorm;
		cplx residualNorm = blas.VectorNorm(r);
		Timers["Error Estimation"].Stop();


		ErrorEstimateList(CurrentArnoldiStep) = residualNorm;
		ErrorEstimate = residualNorm;

		//cout << "Error Estimate = " << ErrorEstimate << " > " << Tolerance << endl;
	}


	MatrixType currentArnoldiMatrix(ArnoldiVectors(blitz::Range(0, CurrentArnoldiStep), blitz::Range::all()));
	VectorType x = HessenbergSolution( blitz::Range(0, CurrentArnoldiStep) );
	//if (HappyBreakdown)
	//{
	//	cout << "Solution = " << x << endl;
	//}
	blas.MultiplyMatrixVector(currentArnoldiMatrix, x, solution);

	Timers["Total"].Stop();
}

template<class T>
void GMRES<T>::PrintStatistics()
{
	if (ProcId == 0)
	{
		cout << "    OpCount      = " << OperatorCount << endl;
		cout << "    InnerProduct = " << InnerProductCount << endl;
		cout << "    ReOrthoCOunt = " << OrthogonalizationCount << endl;
		cout << endl;
		cout << "Timers: " << endl;
		
		for (TimerMap::iterator p=Timers.begin(); p!=Timers.end(); ++p)
		{
			cout << "    " << p->first << " = " << (double)p->second << endl;
		}
	}
}

} //Namespace


#endif


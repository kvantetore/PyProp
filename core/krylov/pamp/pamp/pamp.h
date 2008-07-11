#ifndef PAMP_H
#define PAMP_H

#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include <core/common.h>
#include <core/mpi/mpitraits.h>
#include <core/utility/blitzlapack.h>
#include <core/utility/timer.h>

#include "../../piram/piram/blitzblas.h"
#include "../../piram/piram/functors.h"


namespace pamp
{
using namespace blitz::linalg;

using namespace piram;

/*
 * Main pAMP-class
 */
template <class T>
class pAMP
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

	int PadeOrder;
	int ScalingOrder;

	//Options
	MPI_Comm CommBase;
	bool DisableMPI;

	//Operator class, must implement operator()(VectorType &in, VectorType &out)
	typename OperatorFunctor<T>::Ptr MatrixOperator;
	typename IntegrationFunctor<T, NormType>::Ptr Integration;

	//Constructor
	pAMP() 
	{	
		//Default Params
		PadeOrder = 6;
		ScalingOrder = -1;

		MaxOrthogonalizationCount = 3;
		CommBase = MPI_COMM_WORLD;
		DisableMPI = false;

		Tolerance = sqrt(std::numeric_limits<NormType>::epsilon());
	}

private:
	LAPACK lapack;
	BLAS blas;

	int ProcId;
	int ProcCount;

	//Iteration variables
	/*
	 * Matrices in pAMP are always stored in col-major transposed
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
	MatrixType HessenbergExp;
	MatrixType PadeP, PadeQ, PadeTemp, PadeTemp2;
	IntVectorType PadePivot;
	blitz::Array<double, 1> PadeCoeff;

	VectorType Overlap;
	VectorType Overlap2;
	VectorType Overlap3;
	//Scalars
	int CurrentArnoldiStep;
	bool IsConverged;
	bool HappyBreakdown;
	double Tolerance;
	//...
	
	//Statistics
	int OperatorCount;
	int OrthogonalizationCount;
	int InnerProductCount;
	TimerMap Timers;

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
	void PadeExponential(MatrixType A, MatrixType output, int order);
	void ScalingPadeExponential(MatrixType A, MatrixType output, int scalingOrder, int padeOrder);

	double EstimateMemoryUsage(int matrixSize, int basisSize);
	void Setup();
	void PropagateVector(VectorType input, T dt);
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
		
};

template <class T>
double pAMP<T>::EstimateMemoryUsage(int matrixSize, int basisSize)
{
	double largeSize = matrixSize * (basisSize + 3.0);
	double smallSize = basisSize * (3*basisSize + 7.0);
	return (largeSize + smallSize) * sizeof(T) / (1024.0*1024.0);
}


template <class T>
void pAMP<T>::Setup()
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
	HessenbergMatrix.resize(BasisSize, BasisSize);
	HessenbergExp.resize(BasisSize, BasisSize);
	Overlap.resize(BasisSize);
	Overlap2.resize(BasisSize);
	Overlap3.resize(BasisSize);
	PadeP.resize(BasisSize, BasisSize);
	PadeQ.resize(BasisSize, BasisSize);
	PadeTemp.resize(BasisSize, BasisSize);
	PadeTemp2.resize(BasisSize, BasisSize);
	PadePivot.resize(BasisSize);
	PadeCoeff.resize(30);

	ResetStatistics();

	//Reset arnoldi factorization
	ArnoldiVectors = 0;
	HessenbergMatrix = 0;
	CurrentArnoldiStep = 0;
}

template <class T>
void pAMP<T>::PerformInitialArnoldiStep()
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
}

template <class T>
void pAMP<T>::PerformArnoldiStep()
{
	Timers["Arnoldi Step"].Start();

	//Define some shortcuts
	int j = CurrentArnoldiStep;
	VectorType currentArnoldiVector(ArnoldiVectors(j+1, blitz::Range::all()));
	MatrixType currentArnoldiMatrix(ArnoldiVectors(blitz::Range(0, j+1), blitz::Range::all()));
	VectorType currentOverlap(Overlap(blitz::Range(0, j+1)));
	VectorType currentOverlap2(Overlap2(blitz::Range(0, j+1)));

	//Update the Arnoldi Factorization with the previous Residual
	T beta = CalculateGlobalNorm(Residual);
	if (std::abs(beta) < Tolerance)
	{
		HappyBreakdown = true;
		Timers["Arnoldi Step"].Stop();
		return;
	}
	blas.ScaleVector(Residual, 1.0/beta);
	blas.CopyVector(Residual, currentArnoldiVector);
	HessenbergMatrix(j, j+1) = beta;

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
	//T residualNorm = CalculateGlobalNorm(Residual);
	Timers["Arnoldi Step (Orthogonalization)"].Stop();

	//Update the Hessenberg Matrix
	HessenbergMatrix(j+1, blitz::Range(0, j+1)) = currentOverlap;

	CurrentArnoldiStep++;

	Timers["Arnoldi Step"].Stop();
}

template <class T>
void pAMP<T>::PerformOrthogonalization(T originalNorm, T residualNorm)
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
void pAMP<T>::MultiplyOperator(VectorType &in, VectorType &out)
{
	Timers["Operator"].Start();

	//Update statistics
	OperatorCount++;

	//Perform matrix-vector multiplication
	(*MatrixOperator)(in, out);

	Timers["Operator"].Stop();
}


template<class T>
void pAMP<T>::RecreateArnoldiFactorization()
{
	//Make sure we have a BasisSize-step arnoldi iteration
	while (CurrentArnoldiStep < BasisSize-1 && !HappyBreakdown)
	{
		PerformArnoldiStep();
	}
}


template<class T>
void pAMP<T>::ScalingPadeExponential(MatrixType A, MatrixType output, int padeOrder, int scalingOrder)
{
	Timers["ScalingAndSquaring"].Start();

	if (scalingOrder == -1)
	{
		scalingOrder = ScalingOrder;
	}

	if (scalingOrder < 0)
	{
		//Determine optimal scaling
		double maxValue = blitz::max(abs(A));
		scalingOrder = (int) std::ceil( log2( maxValue ) );
	}
	scalingOrder = std::max(0, scalingOrder);
	//cout << "Using scaling order " << scalingOrder << endl;

	//Scale
	//cout << "Original A = " << A << endl;
	A /= pow(2.0, scalingOrder);
	//cout << "Scaled A = " << A << endl;

	//Calculate Exponential
	Timers["ScalingAndSquaring"].Stop();
	PadeExponential(A, PadeTemp2, padeOrder);
	Timers["ScalingAndSquaring"].Start();

	//Square
	for (int i=0; i<scalingOrder; i++)
	{
		blas.MultiplyMatrixMatrix(PadeTemp2, PadeTemp2, PadeTemp);
		blitz::swap(PadeTemp, PadeTemp2);
	}
	output = PadeTemp2;

	Timers["ScalingAndSquaring"].Stop();
}

template<class T>
void pAMP<T>::PadeExponential(MatrixType A, MatrixType output, int order)
{
	Timers["PadeExponential"].Start();

	if (order == -1)
	{
		order = PadeOrder;
	}
	int N = A.extent(0);

	MatrixType Asqr(PadeTemp2);
	MatrixType P(PadeP);
	MatrixType Q(PadeQ);
	MatrixType temp(PadeTemp);
	IntVectorType pivot(PadePivot);
	blitz::Array<double, 1> coeff(PadeCoeff);

	//Store Asqr = A * A
	blas.MultiplyMatrixMatrix(A, A, Asqr);

	//Calculate coeffs
	coeff(0) = 1;
	for (int i=0; i<order; i++)
	{
		coeff(i+1) = coeff(i) * (double)(order-i) / (double)((i+1) * (2*order-i));
	}

	//Horner eval of P and Q
	P = 0;
	Q = 0;
	for (int i=0; i<N; i++)
	{
		Q(i,i) = coeff(order);
		P(i,i) = coeff(order-1);
	}

	bool odd = true;
	for (int i=order-2; i>=0; i--)
	{
		if (odd)
		{
			blas.MultiplyMatrixMatrix(Q, Asqr, temp);
			blitz::swap(Q, temp);
			for (int row=0; row<N; row++)
			{
				Q(row,row) += coeff(i);
			}
		}
		else
		{
			blas.MultiplyMatrixMatrix(P, Asqr, temp);
			blitz::swap(P, temp);
			for (int row=0; row<N; row++)
			{
				P(row,row) += coeff(i);
			}
		}
		odd = !odd;
	}

	if (odd)
	{
		blas.MultiplyMatrixMatrix(Q, A, temp);
		blitz::swap(Q, temp);
	}
	else
	{
		blas.MultiplyMatrixMatrix(P, A, temp);
		blitz::swap(P, temp);
	}

	Q = Q - P;

	//Do (I + 2*inv(Q)*P)
	int error = lapack.CalculateLUFactorization(Q, pivot);
	if (error != 0)
	{
		cout << "Could solve " << Q << endl;
		cout << "H = " << A << endl;
		throw std::runtime_error("Error in pAMP");
	}
	lapack.SolveGeneralFactored(LAPACK::TransposeNone, Q, pivot, P);
	output = 2.0 * P;
	for (int i=0; i<N; i++)
	{
		output(i,i) += 1;
	}
		
	if (odd)
	{
		output = - output;
	}

	Timers["PadeExponential"].Stop();
}


template<class T>
void pAMP<T>::PropagateVector(VectorType input, T dt)
{
	Timers["Total"].Start();
	Residual = input;

	PerformInitialArnoldiStep();
	RecreateArnoldiFactorization();

	/*
	 * AQQ* = QHQ*
	 *
	 * exp(AQQ*)b = exp(QHQ*)b = Q exp(H) Q* b = Q exp(H) [1, 0, ..., 0]
	 *            = Q exp(H)[0, :]
	 */

	//cout << "H = " << HessenbergMatrix << endl;
	HessenbergMatrix *= - dt * cplx(0., 1.);
	ScalingPadeExponential(HessenbergMatrix, HessenbergExp, -1, -1);

	VectorType q(BasisSize);
	VectorType v(BasisSize);
	q = 0;
	q(0) = 1;

	blas.MultiplyMatrixVector(HessenbergExp, q, v);
	
	blas.MultiplyMatrixVector(ArnoldiVectors, v, input);

	Timers["Total"].Stop();
}

template<class T>
void pAMP<T>::PrintStatistics()
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


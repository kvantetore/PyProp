#include "arpackpropagator.h"

#include <core/wavefunction.h>
#include <core/utility/fortran.h>
#include <core/utility/blitztricks.h>

#include "../krylovcommon.h"
#include "arpack++/pcaupp.h"
#include "arpack++/pceupp.h"
#include "arpack++/caupp.h"
#include "arpack++/ceupp.h"
#include "arpack++/debug.h"

namespace krylov
{

using namespace boost::python;


template<int Rank>
void ArpackPropagator<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	config.Get("krylov_basis_size", BasisSize);
	config.Get("krylov_tolerance", Tolerance);
	config.Get("krylov_eigenvalue_count", EigenvalueCount);
	config.Get("krylov_max_iteration_count", MaxIterationCount);
	config.Get("krylov_use_random_start", RandomStart);
	config.Get("krylov_statistics", PrintStatistics);

	UseParpack = false;
	if (config.HasValue("krylov_use_parpack"))
	{
		config.Get("krylov_use_parpack", UseParpack);
	}
	if (UseParpack)
	{
		cout << "ARPACK: Using Parallel ARPACK" << endl;
	}
}


template<int Rank>
double ArpackPropagator<Rank>::GetMemoryEstimate(const Wavefunction<Rank> &psi)
{
	int matrixSize = psi.GetData().size();
	double size = sizeof(cplx) * matrixSize * (4 + BasisSize);
	return size/(1024*1024);
}



template<int Rank>
void ArpackPropagator<Rank>::Setup(const Wavefunction<Rank> &psi)
{
	int matrixSize = psi.GetData().size();

	//Allocate workspace
	WorkSelect.resize(BasisSize+1);
	WorkVectors.resize(3 * matrixSize);
	WorkVectors2.resize(2 * BasisSize);
	WorkData.resize(3 * sqr(BasisSize) + 5 * BasisSize);
	WorkResidual.resize(matrixSize);
	WorkReal.resize(3 * BasisSize);
	IterationParameters.resize(12);
	IterationParameters = 0;

	Eigenvalues.resize(EigenvalueCount);
	Eigenvectors.resize(BasisSize, matrixSize);
}


template<int Rank>
void ArpackPropagator<Rank>::Solve(object callback, Wavefunction<Rank> &psi, Wavefunction<Rank> &tempPsi)
{
	int matrixSize = psi.GetData().size();
	char* eigenvalueRange = "SR";        //Find eigenvalues of Smallest Magnitude
	bool findEigenvectors = true;         //If we should find eigenvectors or not
	char matrixType = 'I';              //Normal eigenvalueproblem

	int iterationAction  = 0;                      // Action to perform this iteration
	int iterationInfo    = 0;                      // Information from iteration
	blitz::Array<int, 1> iterationPointer(15);     // ?
	int errorNumber      = 0;                      // ierr
	cplx sigma           = 0;                      // ?

	if (RandomStart)
	{
		WorkResidual = 0;
		iterationInfo = 0;
	}
	else
	{
		WorkResidual = MapToRank1(psi.GetData());
		iterationInfo = 1;
	}


	//Set debug info
	int digit = -3;
	int getv0 = 0;
	int aupd = 0;
	int aup2 = 0;
	int aitr = 0;
	int eigt = 0;
	int apps = 0;
	int gets = 0;
	int eupd = 0;
	if (PrintStatistics)
	{
		aupd = 1;
	}
	cTraceOn(digit, getv0, aupd, aup2, aitr, eigt, apps, gets, eupd);


	IterationParameters(1 - 1) = 1;					//shift strategy 
	IterationParameters(3 - 1) = MaxIterationCount; //max number of iterations
	IterationParameters(7 - 1) = 1;                 //mode1 of znaupd is used...

	//Preserve the original psi and tempPsi
	DataArray origPsiData(psi.GetData());
	DataArray origTempPsiData(tempPsi.GetData());

	//---------------------------------------------------------------------------
	//Main loop
	//---------------------------------------------------------------------------
	int i = 0;
	do 
	{
		i = i+1;
		if (UseParpack)
		{
			pcaupp(
				MPI_COMM_WORLD, 
				iterationAction,             // Which action to take after this call (ido)
				matrixType,                  // What type of problem s this (bmat)
				matrixSize,                  // Size of matrix (matrixSize)x(matrixSize) (n)
				eigenvalueRange,             // Which eigenvalues to find (which)
				EigenvalueCount,             // How many eigenvalues to get (nev)
				Tolerance,                   // accuracy on eigenvalues (rol)
				WorkResidual.data(),         // residual (resid)
				BasisSize,                   // number of arnoldi vectors to use (ncv)
				Eigenvectors.data(),         // eigenvectors (v)
				matrixSize,                  // length of eigenvectors (ldv)
				IterationParameters.data(),  // (iparam)
				iterationPointer.data(),     // (ipntr)
				WorkVectors.data(),          // (Workd)
				WorkData.data(),             // (Workl)
				WorkData.size(),             // (lWorkl)
				WorkReal.data(),             // (rWork)
				iterationInfo);
		}
		else
		{
			caupp(
				iterationAction,             // Which action to take after this call (ido)
				matrixType,                  // What type of problem s this (bmat)
				matrixSize,                  // Size of matrix (matrixSize)x(matrixSize) (n)
				eigenvalueRange,             // Which eigenvalues to find (which)
				EigenvalueCount,             // How many eigenvalues to get (nev)
				Tolerance,                   // accuracy on eigenvalues (rol)
				WorkResidual.data(),         // residual (resid)
				BasisSize,                   // number of arnoldi vectors to use (ncv)
				Eigenvectors.data(),         // eigenvectors (v)
				matrixSize,                  // length of eigenvectors (ldv)
				IterationParameters.data(),  // (iparam)
				iterationPointer.data(),     // (ipntr)
				WorkVectors.data(),          // (Workd)
				WorkData.data(),             // (Workl)
				WorkData.size(),             // (lWorkl)
				WorkReal.data(),             // (rWork)
				iterationInfo);
		}
			
		if (abs(iterationAction) == 1)
		{
			//Get the input and output buffer for the matrix mul
			//cout << "Input Index  = " << iterationPointer(0) << endl;
			//cout << "Output Index = " << iterationPointer(1) << endl;
			cplx *inStart = &WorkVectors(iterationPointer(0)-1);
			cplx *outStart = &WorkVectors(iterationPointer(1)-1);

			//Map the data buffers to a a blitz array of correct shape
			DataVector shape = psi.GetData().shape();
			DataVector stride = psi.GetData().stride();
			DataArray inData(inStart, shape, stride, blitz::neverDeleteData);
			DataArray outData(outStart, shape, stride, blitz::neverDeleteData);
			outData = 0;

			//Set psi and tempPsi to point to the data buffers
			psi.SetData(inData);
			tempPsi.SetData(outData);

			//Perform matrix-vector multiplication
			callback(psi, tempPsi);
		}

	} while(abs(iterationAction) == 1);

	//Restore psi and tempPsi
	psi.SetData(origPsiData);
	tempPsi.SetData(origTempPsiData);

	//---------------------------------------------------------------------------
	//Postprocess data   
	//---------------------------------------------------------------------------

	//check for error
	if (iterationInfo < 0) 
	{
		cout << "" << endl;
		cout << "Error during diagonalization, info = " << iterationInfo << endl;
		cout << "" << endl;
	}

	//call zneupd to postprocess data
	//OMG! for en j**** lang liste parametere//
	if (UseParpack)
	{
		pceupp( 
			MPI_COMM_WORLD, 
			findEigenvectors,           // (rvec)
			'A',                        // 
			Eigenvalues.data(),         // (d)
			Eigenvectors.data(),        // (v)
			matrixSize,                 // (ldv)
			sigma,                      // (sigma)
			WorkVectors2.data(),        // (Workev)
			matrixType,                 // (bmat)
			matrixSize,                 // (n)
			eigenvalueRange,            // (which)
			EigenvalueCount,            // (nev)
			Tolerance,                  // (tol)
			WorkResidual.data(),        // (resit)
			BasisSize,                  // (ncv)
			Eigenvectors.data(),        // (v)
			matrixSize,                 // (ldv)
			IterationParameters.data(), // (iparam)
			iterationPointer.data(),    // (ipntr)
			WorkVectors.data(),         // (Workd)
			WorkData.data(),            // (Workl)
			WorkData.size(),            // (lWorkl)
			WorkReal.data(),            // (rWork)
			errorNumber);               // (ierr)
	}
	else
	{
		ceupp(
			findEigenvectors,           // (rvec)
			'A',                        // 
			Eigenvalues.data(),         // (d)
			Eigenvectors.data(),        // (v)
			matrixSize,                 // (ldv)
			sigma,                      // (sigma)
			WorkVectors2.data(),        // (Workev)
			matrixType,                 // (bmat)
			matrixSize,                 // (n)
			eigenvalueRange,            // (which)
			EigenvalueCount,            // (nev)
			Tolerance,                  // (tol)
			WorkResidual.data(),        // (resit)
			BasisSize,                  // (ncv)
			Eigenvectors.data(),        // (v)
			matrixSize,                 // (ldv)
			IterationParameters.data(), // (iparam)
			iterationPointer.data(),    // (ipntr)
			WorkVectors.data(),         // (Workd)
			WorkData.data(),            // (Workl)
			WorkData.size(),            // (lWorkl)
			WorkReal.data(),            // (rWork)
			errorNumber);               // (ierr)
	}	

}

template class ArpackPropagator<1>;
template class ArpackPropagator<2>;
template class ArpackPropagator<3>;
template class ArpackPropagator<4>;

} //Namespace
	

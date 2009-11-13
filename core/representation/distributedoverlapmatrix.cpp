#include "distributedoverlapmatrix.h"
#include "../utility/blitztricks.h"
#include <mpi.h>

#include <Epetra_LinearProblem.h>

#include <Amesos_ConfigDefs.h>
#include <Amesos.h>
#include <Amesos_Superludist.h>

#include <Teuchos_RCP.hpp>

using namespace Teuchos;
using namespace blitz;

template<int Rank>
void DistributedOverlapMatrix<Rank>::SetupRank(Wavefunction<Rank> &srcPsi, int opRank)
{

	//Sanity check: operation rank should be less than rank of wavefunction (and nonzero, duh)
	cout << Rank << " " << opRank << endl;
	assert(opRank < Rank);
	assert(opRank > -1);

	//Assert non-orthogonal rank opRank
	assert (!srcPsi.GetRepresentation()->IsOrthogonalBasis(opRank));

	//Create Epetra map for this rank
	typename Wavefunction<Rank>::Ptr tmpPsi = srcPsi.Copy();
	//Epetra_Map_Ptr wavefunctionMap = CreateWavefunctionMultiVectorEpetraMap<Rank>(tmpPsi, opRank);
	WavefunctionMaps(opRank) = CreateWavefunctionMultiVectorEpetraMap<Rank>(tmpPsi, opRank);

	//Get overlap matrix for opRank (and full col data)
	OverlapMatrix::Ptr overlap = srcPsi.GetRepresentation()->GetGlobalOverlapMatrix(opRank);
	blitz::Array<double, 2> overlapFullCol = overlap->GetOverlapFullCol();

	//Overlap row size
	int numSuperDiagonals = overlap->GetSuperDiagonals();
	int colSizeFull = overlap->GetBasisSize();

	//Epetra CrsMatrix for opRank overlap
	cout << "Allocating Epetra CrsMatrix " << numSuperDiagonals << endl;
	//int numTotalBands = 2 * numSuperDiagonals - 1;
	//Epetra_CrsMatrix_Ptr overlapMatrix = Epetra_CrsMatrix_Ptr( new Epetra_CrsMatrix(Copy, *wavefunctionMap, numTotalBands, true) );

	//Calculate the number of elements per (proc local) row
	int globalStartRow = WavefunctionMaps(opRank)->MinMyGID();
	int globalEndRow = WavefunctionMaps(opRank)->MaxMyGID();
	int numRowsProc = globalEndRow - globalStartRow + 1;
	blitz::Array<int, 1> RowLengths;
	RowLengths.resize(numRowsProc);
	for (int i=globalStartRow; i<=globalEndRow; i++)
	{
		int startCol = std::max(i - numSuperDiagonals, 0);
		int endCol = std::min(i + numSuperDiagonals + 1, colSizeFull);
		RowLengths(i - globalStartRow) = endCol - startCol + 1;
	}

	OverlapMatrices(opRank) = Epetra_CrsMatrix_Ptr( new Epetra_CrsMatrix(Copy, *WavefunctionMaps(opRank), RowLengths.data(), true) );
	OverlapMatrices(opRank)->SetTracebackMode(4);

	/*int procId = 0;
	int procCount = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);
	for (int i=0; i<procCount; i++)
	{
		if (procId == i)
		{
			cout << "Processor " << i << endl;
			cout << "Wavefunction Map " << *WavefunctionMaps(opRank) << endl;
			cout << endl;
		}
	}*/

	//Copy overlap slice corresponding to this proc into CrsMatrix. Since we have the entire overlap
	//matrix (for opRank) on every proc, global and local indices are the same. We only have to
	//determine start and end of row slice for this proc, and the first column index (banded case).
	for (int i=globalStartRow; i<=globalEndRow; i++)
	{
		int startCol = std::max(i - numSuperDiagonals, 0);
		int endCol = std::min(i + numSuperDiagonals + 1, colSizeFull);
		for (int j=startCol; j<endCol; j++)
		{
			int bandIdx = j - i + numSuperDiagonals;

			OverlapMatrices(opRank)->InsertGlobalValues(i, 1, &overlapFullCol(i, bandIdx), &j);
		}
	}

	//Signal end of matrix input
	OverlapMatrices(opRank)->FillComplete();
}


template<int Rank>
shared_ptr<Epetra_MultiVector> DistributedOverlapMatrix<Rank>::SetupMultivector(Wavefunction<Rank> &srcPsi, int opRank)
{
	//Setup Epetra stuff for opRank
	cout << Rank << " " << opRank << endl;
	SetupRank(srcPsi, opRank);

	//Map wavefunction to 3D array (compress before- and after-ranks)
	blitz::Array<cplx, 3> psiData = MapToRank3(srcPsi.GetData(), opRank, 1);
	int beforeSize = psiData.extent(0);
	int opSize = psiData.extent(1);
	int afterSize = psiData.extent(2);
	int otherSize = beforeSize*afterSize;

	//Copy real and imag part of wavefunction into 2D array
	blitz::Array<double, 2> data;
	data.resize(2 * otherSize, opSize);
	for (int i=0; i<beforeSize; i++)
	{
		for (int j=0; j<opSize; j++)
		{
			for (int k=0; k<afterSize; k++)
			{
				data(i*afterSize + k, j) = real( psiData(i,j,k) );
				data(i*afterSize + k + otherSize, j) = imag( psiData(i,j,k) );
			}
		}
	}

	//Number of vectors (in multivector)
	int numVectors = 2 * otherSize; 

	//Create Epetra multivector (view of data)
	Epetra_MultiVector_Ptr srcVec = Epetra_MultiVector_Ptr(new Epetra_MultiVector(Copy, *WavefunctionMaps(opRank), data.data(), opSize, numVectors));

	return srcVec;
}

/*
 * MultiVectorToWavefunction
 *
 */
template<int Rank>
void DistributedOverlapMatrix<Rank>::MultiVectorToWavefunction(Wavefunction<Rank> &psi, Epetra_MultiVector_Ptr vec, int opRank)
{
	//Put result into destPsi
	blitz::Array<cplx, 3> data = MapToRank3(psi.GetData(), opRank, 1);
	int beforeSize = data.extent(0);
	int opSize = data.extent(1);
	int afterSize = data.extent(2);
	int otherSize = beforeSize*afterSize;
	//double *destVecView=0;
	//destVec->ExtractView(&destVecView, opSize);
	//Array<double, 2> data2(destVecView, data.shape(), data.strides(), neverDeleteData)
	blitz::Array<double, 2> data2;
	data2.resize(2 * otherSize, opSize);
	vec->ExtractCopy(data2.data(), opSize);
	for (int i=0; i<beforeSize; i++)
	{
		for (int j=0; j<opSize; j++)
		{
			for (int k=0; k<afterSize; k++)
			{
				double realVal = data2(i*afterSize + k, j);
				double imagVal = data2(i*afterSize + k + otherSize, j);
				data(i,j,k) = realVal + I * imagVal;
			}
		}
	}
}


template<int Rank>
void DistributedOverlapMatrix<Rank>::MultiplyOverlapRank(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &destPsi, int opRank)
{
	SetupRank(srcPsi, opRank);

	//Setup source multivector
	Epetra_MultiVector_Ptr srcVec =  SetupMultivector(srcPsi, opRank);

	//Output multivector (copy of input)
	Epetra_MultiVector_Ptr destVec = Epetra_MultiVector_Ptr(new Epetra_MultiVector(*srcVec));

	//Perform matrix-vector multiplication
	OverlapMatrices(opRank)->Multiply(false, *srcVec, *destVec);

	//Put result multivector into destination wavefunction
	MultiVectorToWavefunction(destPsi, destVec, opRank);
}


template<int Rank>
void DistributedOverlapMatrix<Rank>::SolveOverlapRank(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &destPsi, int opRank)
{
	SetupRank(srcPsi, opRank);
	
	//Setup source multivector
	Epetra_MultiVector_Ptr srcVec =  SetupMultivector(srcPsi, opRank);

	//Output multivector (copy of input)
	Epetra_MultiVector_Ptr destVec = Epetra_MultiVector_Ptr(new Epetra_MultiVector(*srcVec));

	//Set up Epetra LinearProblem with overlap for this rank and input/output multivectors
	Epetra_LinearProblem Problem(OverlapMatrices(opRank).get(), destVec.get(), srcVec.get());

	//Initialize Amesos solver and factory
	Amesos_BaseSolver* Solver;
	Amesos Factory;
	//std::string SolverType = "Amesos_Superludist";
	std::string SolverType = "Amesos_Superludist";
	Solver = Factory.Create(SolverType, Problem);

	//Check that requested solver exists in Amesos
	if (Solver == 0)
	{
		throw std::runtime_error("Specified Amesos solver not available");
	}

	//Setup the parameter list for the solver
	Teuchos::ParameterList List;
	List.set("PrintTiming", true);
	List.set("PrintStatus", true);
	Solver->SetParameters(List);

	//Factorization and solve
	Solver->SymbolicFactorization();
	Solver->NumericFactorization();
	AMESOS_CHK_ERRV(Solver->Solve());

	//Put solution multivector into destination wavefunction
	MultiVectorToWavefunction(destPsi, destVec, opRank);

	Solver->PrintTiming();
}

template class DistributedOverlapMatrix<1>;
template class DistributedOverlapMatrix<2>;
template class DistributedOverlapMatrix<3>;
template class DistributedOverlapMatrix<4>;


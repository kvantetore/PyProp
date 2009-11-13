#include "distributedoverlapmatrix.h"
#include "../trilinos/pyprop_epetra.h"
#include <Epetra_MultiVector.h>
#include "../utility/blitztricks.h"

#include <mpi.h>

typedef shared_ptr<Epetra_MultiVector> Epetra_MultiVector_Ptr;

template<int Rank>
void DistributedOverlapMatrix<Rank>::MultiplyOverlapRank(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &destPsi, int opRank)
{
	using namespace blitz;

	//Sanity check: operation rank should be less than rank of wavefunction (and nonzero, duh)
	assert(opRank < Rank);
	assert(opRank > -1);

	//Assert non-orthogonal rank opRank
	assert (!srcPsi.GetRepresentation()->IsOrthogonalBasis(opRank));

	//Create Epetra map for this rank
	typename Wavefunction<Rank>::Ptr tmpPsi = srcPsi.Copy();
	Epetra_Map_Ptr wavefunctionMap = CreateWavefunctionMultiVectorEpetraMap<Rank>(tmpPsi, opRank);

	//Map wavefunction to 3D array (compress before- and after-ranks)
	Array<cplx, 3> psiData = MapToRank3(srcPsi.GetData(), opRank, 1);
	int beforeSize = psiData.extent(0);
	int opSize = psiData.extent(1);
	int afterSize = psiData.extent(2);
	int otherSize = beforeSize*afterSize;

	//Copy real and imag part of wavefunction into 2D array
	Array<double, 2> data;
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
	Epetra_MultiVector_Ptr srcVec = Epetra_MultiVector_Ptr(new Epetra_MultiVector(Copy, *wavefunctionMap, data.data(), opSize, numVectors));

	//Output multivector (copy of input)
	Epetra_MultiVector_Ptr destVec = Epetra_MultiVector_Ptr(new Epetra_MultiVector(*srcVec));

	/*
	 *
	 * Setup Epetra overlap matrix
	 *
	 */

	//Get overlap matrix for opRank (and full col data)
	OverlapMatrix::Ptr overlap = srcPsi.GetRepresentation()->GetGlobalOverlapMatrix(opRank);
	Array<double, 2> overlapFullCol = overlap->GetOverlapFullCol();

	//Overlap row size
	int numSuperDiagonals = overlap->GetSuperDiagonals();
	int numTotalBands = 2 * numSuperDiagonals - 1;
	int fullRowSize = overlap->GetBasisSize();

	//Epetra CrsMatrix for opRank overlap
	cout << "Allocating Epetra CrsMatrix " << numSuperDiagonals << endl;
	Epetra_CrsMatrix_Ptr overlapMatrix = Epetra_CrsMatrix_Ptr( new Epetra_CrsMatrix(Copy, *wavefunctionMap, numTotalBands, true) );
	overlapMatrix->SetTracebackMode(4);

	//Copy overlap slice corresponding to this proc into CrsMatrix. Since we have the entire overlap
	//matrix (for opRank) on every proc, global and local indices are the same. We only have to
	//determine start and end of row slice for this proc, and the first column index (banded case).

	int procId = 0;
	int procCount = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);
	for (int i=0; i<procCount; i++)
	{
		if (procId == i)
		{
			cout << "Processor " << i << endl;
			cout << "Wavefunction Map " << *wavefunctionMap << endl;
			cout << endl;
		}
	}


	int globalStartRow = wavefunctionMap->MinMyGID();
	int globalEndRow = wavefunctionMap->MaxMyGID();
	for (int i=globalStartRow; i<=globalEndRow; i++)
	{
		int startCol = std::max(i - numSuperDiagonals, 0);
		int endCol = std::min(i + numSuperDiagonals + 1, fullRowSize);
		for (int j=startCol; j<endCol; j++)
		{
			int bandIdx = j - i + numSuperDiagonals;

			overlapMatrix->InsertGlobalValues(i, 1, &overlapFullCol(i, bandIdx), &j);
		}
	}

	//Signal end of matrix input
	overlapMatrix->FillComplete();

	//Perform matrix-vector multiplication
	overlapMatrix->Multiply(false, *srcVec, *destVec);

	//Put result into destPsi
	Array<cplx, 3> destPsiData = MapToRank3(destPsi.GetData(), opRank, 1);
	//double *destVecView=0;
	//destVec->ExtractView(&destVecView, opSize);
	//Array<double, 2> data2(destVecView, data.shape(), data.strides(), neverDeleteData)
	Array<double, 2> data2;
	data2.resize(2 * otherSize, opSize);
	destVec->ExtractCopy(data2.data(), opSize);
	for (int i=0; i<beforeSize; i++)
	{
		for (int j=0; j<opSize; j++)
		{
			for (int k=0; k<afterSize; k++)
			{
				double realVal = data2(i*afterSize + k, j);
				double imagVal = data2(i*afterSize + k + otherSize, j);
				destPsiData(i,j,k) = realVal + I * imagVal;
			}
		}
	}
}


template class DistributedOverlapMatrix<1>;
template class DistributedOverlapMatrix<2>;
template class DistributedOverlapMatrix<3>;
template class DistributedOverlapMatrix<4>;

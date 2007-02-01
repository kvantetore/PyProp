
#include "distributedmodel.h"

#ifndef SINGLEPROC
#include <mpi.h>
#endif

using namespace blitz;

template<int Rank>
DistributedModel<Rank>::DistributedModel() :
	MpiTime(0.0),
	TotalTime(0.0),
	BufferName1(0),
	BufferName2(0)
{
	
#ifndef SINGLEPROC
	MPI_Comm_rank(MPI_COMM_WORLD, &this->ProcId);
	MPI_Comm_size(MPI_COMM_WORLD, &this->ProcCount);
#else
	this->ProcId = 0;
	this->ProcCount = 1;
#endif
	
	std::cout << "Proc " << ProcId << "/" << ProcCount << std::endl;
}

/**
Change the distributed representation according to Model
*/
template<int Rank>
void DistributedModel<Rank>::ChangeRepresentation(Wavefunction<Rank> &psi)
{
	TotalTimer.start();

	int srcName;
	int dstName;
	if (BufferName1 == -1)
	{
		BufferName1 = psi.GetActiveBufferName();
		if (ProcCount > 1)
		{
			BufferName2 = psi.AllocateData(	psi.GetData().shape() );
		}
		else
		{
			BufferName2 = BufferName1;
		}
	}

	if (psi.GetActiveBufferName() == BufferName1) 
	{
		srcName = BufferName1;
		dstName = BufferName2;
	}
	else if (psi.GetActiveBufferName() == BufferName2)
	{
		srcName = BufferName2;
		dstName = BufferName1;
	}
	else
	{
		throw std::runtime_error("DistributedModel::ChangeRepresentation: Active buffer has changed!");
	}
	
	blitz::Array<cplx, Rank> src = psi.GetData(srcName);
	blitz::Array<cplx, Rank> dst = psi.GetData(dstName);

	//Choose direction
	TransposeDirection direction = TRANSPOSE_FORWARD;
	if (
		(Model == TRANSPOSE_SEMI && GetDistributedRank(psi) != 0) 
		|| (Model == TRANSPOSE_SEMI && GetDistributedRank(psi) != 0)
	)
	{
		direction = TRANSPOSE_BACKWARD;
	}
	
	//Transpose into tempdata
	if (Model == TRANSPOSE_SEMI)
	{
		psi.DistributedRank = TransposeSemi(src, dst, psi.DistributedRank, direction);
	}
	else
	{
		Transpose(src, dst, direction);
		psi.DistributedRank = dst.ordering(Rank-1);
	}

	//Make the transposed data active
	psi.SetActiveBuffer(dstName);
	
	TotalTimer.stop();
	TotalTime += (double)TotalTimer;
}

template<int Rank>
int DistributedModel<Rank>::TransposeSemi(blitz::Array<cplx,Rank> &A, blitz::Array<cplx,Rank> &B, int distributedRank, TransposeDirection direction)
{
	const blitz::TinyVector<int, Rank> rankOrder = A.ordering();
	const blitz::TinyVector<int, Rank> extent = A.shape();

	/* 
	Create a 2D view of the entire matrix, where the dimensions with stride < MaxStride are
	treated as one huge dimension. In other words, the 2D matrix has a huge number of coloumns 
	(N = N_0 * N_1 * ... * N_(N-2)) and not as many rows (M = N_(N-1)).
	
	When we transpose this 2D matrix, we end up dividing the dimension of the second 
	highest stride when we're going forward, and the highest stride when we are going
	backward.
	*/	
	int sizeN = 1;
	int sizeM = 1;
	
	int curDistributedRank = -1;
	int newDistributedRank = -1;
	
	if (direction == TRANSPOSE_FORWARD) 
	{
		curDistributedRank = rankOrder(Rank - 1);
		newDistributedRank = rankOrder(Rank - 2);
	}
	else
	{
		curDistributedRank = rankOrder(Rank - 2);
		newDistributedRank = rankOrder(Rank - 1);
	}
	if (curDistributedRank != distributedRank) 
	{
		std::cout << "transposeSemi(): Could not transpose rank " << distributedRank 
			  << ". It is not the correct rank to transpose at this time." << std::endl;
		throw runtime_error("Error during transposeSemi()");
	}
	
	sizeM = extent(curDistributedRank);
	sizeN = 1;
	for (int curRank=0; curRank<Rank; curRank++)
	{
		if (curRank != curDistributedRank)
		{
			sizeN *= extent(curRank);
		}
	}
	
	/* 
	A-matrix is in column-first (C style) ordering.
	*/
	blitz::TinyVector<int,2> shape2dA(sizeM, sizeN);
	blitz::TinyVector<int,2> stride2dA(sizeN, 1);
	
	/*
	In this version (transposeSemi), we do not do a full transpose, 
	but leave the bMatrix in C style ordering. That way we save a LOT of 
	cache misses during transpose2d
	*/
	blitz::TinyVector<int,2> shape2dB(sizeM * ProcCount, sizeN / ProcCount);
	blitz::TinyVector<int,2> stride2dB(shape2dB(secondDim), 1);

	if (direction == TRANSPOSE_BACKWARD)
	{
		blitz::TinyVector<int, 2> temp;
		temp = shape2dA;
		shape2dA = shape2dB;
		shape2dB = temp;
		
		temp = stride2dA;
		stride2dA = stride2dB;
		stride2dB = temp;		
	}


	blitz::Array<cplx, 2> aMatrix(A.data(), shape2dA, stride2dA, neverDeleteData);
	blitz::Array<cplx, 2> bMatrix(B.data(), shape2dB, stride2dB, neverDeleteData);

/*	
	for (int i=0;i<ProcCount;i++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (ProcId == i) 
		{
			cout << "aMatrix (" << ProcId << ") " 
			     << aMatrix.shape() << endl 
			     << aMatrix.stride() << endl;
	   		     // << aMatrix << endl
			     ;
		}
	}
*/
	//Perform a 2D transpose A into B
	if (direction == TRANSPOSE_FORWARD)
	{
		Transpose2d(aMatrix, bMatrix);
	}
	else
	{
		Transpose2dBack(aMatrix, bMatrix);
	}

/*	
	for (int i=0;i<ProcCount;i++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (ProcId == i) 
		{
			cout << "bMatrix (" << ProcId << ") " 
			     << bMatrix.shape() << endl 
			     << bMatrix.stride() << endl 
			     //<< bMatrix << endl
			     ;
		}
	}
*/

	//Make sure B is properly set up.
	blitz::TinyVector<int,Rank> bShape;
	
	//Set up the shape of B
	bShape = extent;
	//The previously distributed rank was divided by ProcCount
	bShape(curDistributedRank) *= ProcCount;
	//The currently distributed rank is now divided by ProcCount
	bShape(newDistributedRank) /= ProcCount;
	
	//Set up B
	if (product(bShape) != product(extent)) 
	{
		cout << "Invalid b-shape " << bShape << ", " << extent << endl;
		exit(-1);
	}
	B.changeShape(bShape);
	
	return newDistributedRank;
}


template<int Rank>
void DistributedModel<Rank>::Transpose(blitz::Array<cplx,Rank> &A, blitz::Array<cplx,Rank> &B, TransposeDirection direction)
{
	const blitz::TinyVector<int, Rank> rankOrder = A.ordering();
	const blitz::TinyVector<int, Rank> extent = A.shape();

	/* 
	Create a 2D view of the entire matrix, where the dimensions with stride < MaxStride are
	treated as one huge dimension. In other words, the 2D matrix has a huge number of coloumns 
	(N = N_0 * N_1 * ... * N_(N-2)) and not as many rows (M = N_(N-1)).
	
	When we transpose this 2D matrix, we end up dividing the dimension of the second 
	highest stride
	*/	
	int sizeN = 1;
	int sizeM = 1;
	
	if (direction == TRANSPOSE_FORWARD) 
	{
		//firstDim (M) has the extent of the Rank with the largest stride
		sizeM = extent(rankOrder(Rank-1));
		
		//secondDim (N) has the extent of all other ranks combined 
		sizeN = 1;
		for (int curRank=0; curRank<Rank-1; curRank++)
		{
			sizeN *= extent(rankOrder(curRank));
		}
	}
	else
	{
		//firstDim (M) has the extent of the Rank with the largest stride
		sizeN = extent(rankOrder(0));
		
		//secondDim (N) has the extent of all other ranks combined 
		sizeM = 1;
		for (int curRank=1; curRank<Rank; curRank++)
		{
			sizeM *= extent(rankOrder(curRank));
		}
	}
	
	/* 
	A-matrix is in column-first (C style) ordering.
	*/
	blitz::TinyVector<int,2> shape2dA(sizeM, sizeN);
	blitz::TinyVector<int,2> stride2dA(sizeN, 1);
	blitz::Array<cplx, 2> aMatrix(A.data(), shape2dA, stride2dA, neverDeleteData);
	
	/*
	B-matrix is in row-first (fortran style) ordering.
	This way, we get a full transpose, and end up with the
	Previously MaxStrided rank (the distributed rank in the A-representation)
	getting stride = 1. 
	When direction is TRANSPOSE_NONE, we shall not change any strides, and therefore 
	leave B in C-style ordering.
	*/
	blitz::TinyVector<int,2> shape2dB(sizeM * ProcCount, sizeN / ProcCount);
	blitz::TinyVector<int,2> stride2dB;
	stride2dB = 1, shape2dB(firstDim);
	blitz::Array<cplx, 2> bMatrix(B.data(), shape2dB, stride2dB, neverDeleteData);
	
	//Perform a 2D transpose A into B
	Transpose2d(aMatrix, bMatrix);
	
	//Make sure B is properly set up.
	blitz::TinyVector<int,Rank> bShape;
	blitz::TinyVector<int,Rank> bRankOrder;
	
	/*
	Set up the ordering of B
	
	Transposing A to B has the effect of making the highest strided rank in A,
	the lowest strided rank in B, and increase the order of stride by one on all
	other ranks.  
	i.e. order = thirdDim, secondDim, firstDim 
	-> order = firstDim, thirdDim, secondDim
	(or the other way around, if direction is TRANSPOSE_BACKWARD)
	*/
	for (int curRank=0; curRank<Rank; curRank++)
	{
		bRankOrder(curRank) = rankOrder( (curRank + Rank + direction ) % Rank );
	}
	
	//Set up the shape of B
	bShape = extent;
	//The previously highestly strided rank was divided by ProcCount
	bShape(rankOrder(Rank - 1)) *= ProcCount;
	//The currently highestly strided rank is going to be divided by ProcCount
	bShape(bRankOrder(Rank - 1)) /= ProcCount;
		
	//Set up B
	B.changeOrdering(bRankOrder);
	if (product(bShape) != product(extent)) 
	{
		std::cout << "Invalid b-shape " << bShape << ", " << extent << endl;
		exit(-1);
	}
	B.changeShape(bShape);
}

template<class T>
void blitz_copy(blitz::Array<T,2> &src, int srcStartM, int srcStartN, blitz::Array<T,2> &dst, int dstStartM, int dstStartN, int countM, int countN)
{
	if (src.stride(0) < src.stride(1) || dst.stride(0) < dst.stride(1))
	{
		cout << "Can only perform fast array copy when both arrays are C-style ordered" << endl;
		throw runtime_error("distributedmodel.cpp, blitz_copy() Invalide strides");
	}
	
	size_t rowsize = countN * sizeof(T);

	T* srcPtr = (T*)src.data();
	srcPtr += srcStartM * src.stride(0) + srcStartN * src.stride(1);
			
	T* dstPtr = (T*)dst.data();
	dstPtr += dstStartM * dst.stride(0) + dstStartM * dst.stride(1);
	
	for (int m=0; m<countM; m++)
	{
		memcpy(dstPtr, srcPtr, rowsize);
		dstPtr += dst.stride(0);
		srcPtr += src.stride(0);
	}	
}


template<int Rank>
void DistributedModel<Rank>::Transpose2d(blitz::Array<cplx,2> &A, blitz::Array<cplx,2> &B)
{
	
	//Common variables
	int procSizeM = A.extent(firstDim);
	int procSizeN = A.extent(secondDim) / ProcCount;
#ifndef SINGLEPROC
	int blockSize = procSizeM * procSizeN * sizeof(cplx) / sizeof(double);
	bool useRecieveBuffer = false;
	
	Array<cplx, 2> tempA;
	Array<cplx, 2> tempB;

	//If we have only one proc, we really don't need the buffers
	if (ProcCount > 0) 
	{
		//Resize tempA-buffer
		blitz::TinyVector<int, 2> newTempShapeA(procSizeM, procSizeN);
		if (blitz::product(newTempShapeA) == tempB.size())
		{
			tempA.changeShape(newTempShapeA);
		}
		else
		{
			tempA.resize(newTempShapeA);
		}

		//If B is not in C-ordering, we must use a temporary recieve buffer
		if (B.stride(0) < B.stride(1))
		{
			useRecieveBuffer = true;

			//Resize tempB-buffer
			blitz::TinyVector<int, 2> newTempShapeB(B.extent(firstDim) / ProcCount, B.extent(secondDim));
			if (blitz::product(newTempShapeB) == tempB.size())
			{
				tempB.changeShape(newTempShapeB);
			}
			else
			{
				tempB.resize(newTempShapeB);
			}
		}
	}
#endif	

	//Copy diagonal elements
	//Find the range of M in B which is on the diagonal
	//NB! blitz::Range has the non-intuitive behavior of including
	//The last index in the range, hence the "-1"
	int bDiagonalStartM = ProcId * procSizeM;
	int bDiagonalEndM = bDiagonalStartM + procSizeM;
	Range bDiagonalRangeM(bDiagonalStartM, bDiagonalEndM-1);
	
	//Find the range of N in A which is on the diagonal
	int aDiagonalStartN = ProcId * procSizeN;
	int aDiagonalEndN = aDiagonalStartN + procSizeN;
	Range aDiagonalRangeN(aDiagonalStartN, aDiagonalEndN-1);
	
	//Copy and transpose diagonal elements
	//(B has a different storage order than A, which leads to a transpose)
	B(bDiagonalRangeM, Range::all()) = A(Range::all(), aDiagonalRangeN); 
//	blitz_copy(A, 0, aDiagonalStartN, B, bDiagonalStartM, 0, procSizeM, procSizeN);
	
	/* 
	blocks are numbered per proc as starting on 0 on the diagonal,
	and cycicly to the left (for A) and down (for B).
	For a 4x4 block grid we have
	A: Proc0: (0) (1) (2) (3)
	Proc1: (3) (0) (1) (2)
	Proc2: (2) (3) (0) (1)
	Proc3: (1) (2) (3) (0)

	B:         p0  p1  p2  p3
		(0) (3) (2) (1)
		(1) (0) (3) (2)
		(2) (1) (0) (3)
		(3) (2) (1) (0)
	*/

	//We have already copied the diagonal (block 0), and have ProcCount-1 to go
#ifndef SINGLEPROC
	for (int aBlock=1; aBlock<ProcCount; aBlock++)
	{
		// aBlock and bBlock are at the same distance from the 
		// diagonal, only in the opposite directions
		int bBlock = (ProcCount - aBlock) % ProcCount;
		
		// We need to recieve from the proc that will try to send to us
		// during this value of aBlock, otherwise we will end up with a
		// deadlock!
		int toProc = (aBlock + ProcId) % ProcCount;
		int fromProc = (bBlock + ProcId) % ProcCount;
		
		int aGridIndex = toProc;
		int bGridIndex = fromProc;
		
		//Construct the array range in N for A
		blitz::Range aBlockRangeN(aGridIndex * procSizeN, (aGridIndex + 1) * procSizeN - 1);
		blitz::Range bBlockRangeM(bGridIndex * procSizeM, (bGridIndex + 1) * procSizeM - 1);
		
		//Copy data from A into a temporary array, in order to make the data continous
		//This could perhaps be replaced by using MPI-packing
		
		tempA = A(Range::all(), aBlockRangeN);
//		blitz_copy(A, 0, aGridIndex * procSizeN, tempA, 0, 0, procSizeM, procSizeN);

		
		//We recieve the data in the same order as we intend to store it in B, and 
		//can therefore make MPI store it there directly.
		double* recieveBuffer;
		if (useRecieveBuffer) {
			recieveBuffer = reinterpret_cast<double*>(tempB.data());
		}
		else
		{
			recieveBuffer = reinterpret_cast<double*>(&B(bGridIndex * procSizeM, 0));
		}

		MpiTimer.start();
		MPI_Request request;
		MPI_Irecv(recieveBuffer, blockSize, MPI_DOUBLE, fromProc, 0, MPI_COMM_WORLD, &request);
		MPI_Send(tempA.data(), blockSize, MPI_DOUBLE, toProc, 0, MPI_COMM_WORLD);
		MPI_Status status;
		MPI_Wait(&request, &status);
		MpiTimer.stop();
		MpiTime += (double) MpiTimer;
		
		if (useRecieveBuffer)
		{
			//Copy and transpose data into B
			B(bBlockRangeM, Range::all()) = tempB;
		}
	}
#endif //SINGLEPROC
}


template<int Rank>
void DistributedModel<Rank>::Transpose2dBack(blitz::Array<cplx,2> &A, blitz::Array<cplx,2> &B)
{
	//Common variables
	int procSizeM = A.extent(firstDim) / ProcCount;
	int procSizeN = A.extent(secondDim);
#ifndef SINGLEPROC
	int blockSize = procSizeM * procSizeN * sizeof(cplx) / sizeof(double);
	//bool useRecieveBuffer = false;

	//If we have only one proc, we really don't need the buffers
	if (ProcCount > 1) 
	{
		blitz::TinyVector<int, 2> newTempShape(B.extent(firstDim), B.extent(secondDim) / ProcCount);
		if (blitz::product(newTempShape) == tempB.size())
		{
			tempB.changeShape(newTempShape);
		}
		else
		{
			tempB.resize(newTempShape);
		}
	}
#endif 
	
	//Copy diagonal elements
	//Find the range of N in B which is on the diagonal
	//NB! blitz::Range has the non-intuitive behavior of including
	//The last index in the range, hence the "-1"
	int bDiagonalStartN = ProcId * procSizeN;
	int bDiagonalEndN = bDiagonalStartN + procSizeN;
	Range bDiagonalRangeN(bDiagonalStartN, bDiagonalEndN-1);
	
	//Find the range of M in A which is on the diagonal
	int aDiagonalStartM = ProcId * procSizeM;
	int aDiagonalEndM = aDiagonalStartM + procSizeM;
	Range aDiagonalRangeM(aDiagonalStartM, aDiagonalEndM-1);
	
	//Copy without transpose diagonal elements
	B(Range::all(), bDiagonalRangeN) = A(aDiagonalRangeM, Range::all()); 
//	blitz_copy(A, aDiagonalStartM, 0, B, 0, bDiagonalStartN, procSizeM, procSizeN);
	
	/* 
	blocks are numbered per proc as starting on 0 on the diagonal,
	and cycicly to the left (for B) and down (for A).
	For a 4x4 block grid we have
	B: Proc0: (0) (1) (2) (3)
	   Proc1: (3) (0) (1) (2)
	   Proc2: (2) (3) (0) (1)
	   Proc3: (1) (2) (3) (0)

	A:         p0  p1  p2  p3
	          (0) (3) (2) (1)
	          (1) (0) (3) (2)
	          (2) (1) (0) (3)
	          (3) (2) (1) (0)
	*/

#ifndef SINGLEPROC
	//We have already copied the diagonal (block 0), and have ProcCount-1 to go
	for (int aBlock=1; aBlock<ProcCount; aBlock++)
	{
		// aBlock and bBlock are at the same distance from the 
		// diagonal, only in the opposite directions
		int bBlock = (ProcCount - aBlock) % ProcCount;
		
		// We need to recieve from the proc that will try to send to us
		// during this value of aBlock, otherwise we will end up with a
		// deadlock!
		int toProc = (aBlock + ProcId) % ProcCount;
		int fromProc = (bBlock + ProcId) % ProcCount;
		
		int aGridIndex = toProc;
		int bGridIndex = fromProc;
		
		//Construct the array range in N for A
		Range bBlockRangeN(bGridIndex * procSizeN, (bGridIndex + 1) * procSizeN - 1);
			
		//when we're going back, A is contigous, and B needs upacking
		double* sendBuffer = reinterpret_cast<double*>(&A(aGridIndex * procSizeM, 0));
		double* recieveBuffer = reinterpret_cast<double*>(tempB.data());
	
		MpiTimer.start();
		MPI_Request request;
		MPI_Irecv(recieveBuffer, blockSize, MPI_DOUBLE, fromProc, 0, MPI_COMM_WORLD, &request);
		MPI_Send(sendBuffer, blockSize, MPI_DOUBLE, toProc, 0, MPI_COMM_WORLD);
		MPI_Status status;
		MPI_Wait(&request, &status);
		MpiTimer.stop();
		MpiTime += (double) MpiTimer;
		
		B(Range::all(), bBlockRangeN) = tempB;
//		blitz_copy(tempB, 0, 0, B, 0, bGridIndex*procSizeN, procSizeM, procSizeN);
	}
#endif

}


template<int Rank>
void DistributedModel<Rank>::ApplyConfigSection(const ConfigSection &cfg)
{
	int model = -1;
	cfg.Get("transpose_model", model);
	
	this->Model = (TransposeModel)model;
	std::cout << "Using transpose model " << this->Model << std::endl;
}

template class DistributedModel<1>;
template class DistributedModel<2>;
template class DistributedModel<3>;
template class DistributedModel<4>;
template class DistributedModel<5>;
template class DistributedModel<6>;
template class DistributedModel<7>;
template class DistributedModel<8>;
template class DistributedModel<9>;


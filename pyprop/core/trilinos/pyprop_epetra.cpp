#ifdef PYPROP_USE_TRILINOS
#include "pyprop_epetra.h"

/*
 * Creates an Epetra_Comm object from a DistributedModel. It should
 * use the distributed model to get an MPI_Communicator, but right now
 * DistributedModel uses MPI_COMM_WORLd.
 */
template<int Rank>
Epetra_Comm_Ptr CreateDistributedModelEpetraComm(typename DistributedModel<Rank>::Ptr distr)
{
#ifdef EPETRA_MPI
	//should get MPI_COMM from distributed model
	return shared_ptr<Epetra_MpiComm>( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
	assert(distr->IsSingleProc());
	return shared_ptr<Epetra_SerialComm>( new Epetra_SerialComm() );
#endif
}

template Epetra_Comm_Ptr CreateDistributedModelEpetraComm<1>(DistributedModel<1>::Ptr distr);
template Epetra_Comm_Ptr CreateDistributedModelEpetraComm<2>(DistributedModel<2>::Ptr distr);
template Epetra_Comm_Ptr CreateDistributedModelEpetraComm<3>(DistributedModel<3>::Ptr distr);
template Epetra_Comm_Ptr CreateDistributedModelEpetraComm<4>(DistributedModel<4>::Ptr distr);


/*
 * Creates an Epetra_Comm object from a DistributedModel.
 *
 *   distr: Distributed model in use (possible multi-dim proc grid)
 *   rank: which rank in the wavefunction to use when creating Epetra comm. If this rank is
 *         not distributed, we return a SerialComm for now.
 *
 */
template<int Rank>
Epetra_Comm_Ptr CreateDistributedModelEpetraComm(typename DistributedModel<Rank>::Ptr distr, int wavefunctionRank)
{

#ifdef EPETRA_MPI

	typedef blitz::Array<MPI_Comm, 1> ProcVectorComm;

	//Check if wavefunctionRank is distributed; otherwise no communication is needed
	if (distr->IsDistributedRank(wavefunctionRank))
	{
		//get MPI_COMM from distributed model, for given proc grid rank
		MPI_Comm groupCommRank = distr->GetGroupCommRank(wavefunctionRank);
		return shared_ptr<Epetra_MpiComm>( new Epetra_MpiComm(groupCommRank) );
	}
	else
	{
		return shared_ptr<Epetra_SerialComm>( new Epetra_SerialComm() );
	}
	

#else
	assert(distr->IsSingleProc());
	return shared_ptr<Epetra_SerialComm>( new Epetra_SerialComm() );

#endif
}

template Epetra_Comm_Ptr CreateDistributedModelEpetraComm<1>(DistributedModel<1>::Ptr distr, int procRank);
template Epetra_Comm_Ptr CreateDistributedModelEpetraComm<2>(DistributedModel<2>::Ptr distr, int procRank);
template Epetra_Comm_Ptr CreateDistributedModelEpetraComm<3>(DistributedModel<3>::Ptr distr, int procRank);
template Epetra_Comm_Ptr CreateDistributedModelEpetraComm<4>(DistributedModel<4>::Ptr distr, int procRank);


/*
 * Creates an Epetra_Map corresponding to how a DistributedModel distributes a wavefunction, 
 * for a single rank (opRank). That is, we group together the procs that has the remaining 
 * parts of the wavefunction data for given opRank. Thus, only local <-> global index map
 * for this rank needs be constructed.
 *
 */
template<int Rank>
Epetra_Map_Ptr CreateWavefunctionMultiVectorEpetraMap(typename Wavefunction<Rank>::Ptr psi, int opRank)
{
	typedef blitz::Array<cplx, Rank> ArrayType;
	
	//Create a comm object from wavefunction distribution. Only involves procs with
	//data for current operation rank
	typename DistributedModel<Rank>::Ptr distr = psi->GetRepresentation()->GetDistributedModel();
	Epetra_Comm_Ptr comm = CreateDistributedModelEpetraComm<Rank>(distr, opRank);

	//Local size of wavefunction (opRank)
	ArrayType data = psi->GetData();
	int localSize = data.extent(opRank);

	//Get global start idx for opRank
	int rankSizeGlobal = psi->GetRepresentation()->GetFullShape()(opRank);
	int globalStartIdx = distr->GetLocalStartIndex(rankSizeGlobal, opRank);

	//Create index map for local elements (local idx -> global idx) 
	int* myElems = new int[localSize];

	//Setup local->global index map
	for (int i=0; i<localSize; i++)	
	{
		myElems[i] = i + globalStartIdx;
	}

	return Epetra_Map_Ptr( new Epetra_Map(-1, localSize, myElems, 0, *comm) );
}
template Epetra_Map_Ptr CreateWavefunctionMultiVectorEpetraMap<1>(Wavefunction<1>::Ptr psi, int opRank);
template Epetra_Map_Ptr CreateWavefunctionMultiVectorEpetraMap<2>(Wavefunction<2>::Ptr psi, int opRank);
template Epetra_Map_Ptr CreateWavefunctionMultiVectorEpetraMap<3>(Wavefunction<3>::Ptr psi, int opRank);
template Epetra_Map_Ptr CreateWavefunctionMultiVectorEpetraMap<4>(Wavefunction<4>::Ptr psi, int opRank);


/*
 * Creates an Epetra_Map corresponding to how a DistributedModel
 * distributes a wavefunction.
 *
 * This map can be used to create Epetra_Vectors that are array-compatible
 * with Wavefunctions, i.e. has the same memory layout and distribution
 */
template<int Rank>
Epetra_Map_Ptr CreateWavefunctionEpetraMap(typename Wavefunction<Rank>::Ptr psi)
{
	typedef blitz::Array<cplx, Rank> ArrayType;
	

	ArrayType data = psi->GetData();
	typename DistributedModel<Rank>::Ptr distr = psi->GetRepresentation()->GetDistributedModel();
	Epetra_Comm_Ptr comm = CreateDistributedModelEpetraComm<Rank>(distr);

	int localSize = data.size();
	int* myElems = new int[localSize*2];

	//Get strides of the global Rank-dimensional array
	blitz::TinyVector<int, Rank> globalShape = psi->GetRepresentation()->GetFullShape();
	blitz::TinyVector<int, Rank> globalStrides;
	globalStrides(Rank-1) = 1;
	for (int rank=Rank-2; rank>=0; rank--)
	{
		globalStrides(rank) = globalStrides(rank+1)	* globalShape(rank+1);
	}

	//Get maps to global indices for each rank
	blitz::TinyVector<blitz::Range, Rank> globalIndices;
	for (int rank=0; rank<Rank; rank++)
	{
		globalIndices(rank) = distr->GetLocalIndexRange(globalShape(rank), rank);
	}
	
	typename blitz::Array<cplx, Rank>::iterator it = data.begin();
	for (int linearCount=0; linearCount<data.size(); linearCount++)
	{
		//Get global position for the current element
		int globalPos = 0;
		for (int rank=0; rank<Rank; rank++)
		{
			int rankPos = it.position()(rank);
			globalPos += globalStrides(rank) * globalIndices(rank)(rankPos);
		}

		//real component
		myElems[2*linearCount] = 2*globalPos;
		//imaginary component
		myElems[2*linearCount+1] = 2*globalPos + 1;
			
		++it;	
	}

	return Epetra_Map_Ptr( new Epetra_Map(-1, localSize*2, myElems, 0, *comm) );
}
template Epetra_Map_Ptr CreateWavefunctionEpetraMap<1>(Wavefunction<1>::Ptr psi);
template Epetra_Map_Ptr CreateWavefunctionEpetraMap<2>(Wavefunction<2>::Ptr psi);
template Epetra_Map_Ptr CreateWavefunctionEpetraMap<3>(Wavefunction<3>::Ptr psi);
template Epetra_Map_Ptr CreateWavefunctionEpetraMap<4>(Wavefunction<4>::Ptr psi);

/*
 * Creates an Epetra_Vector that is a view into a wavefunction data.
 *   - the wavefunctionMap should be created with CreateWavefunctionEpetraMap(...)
 */
template<int Rank>
Epetra_Vector_Ptr CreateWavefunctionEpetraView(typename Wavefunction<Rank>::Ptr psi, Epetra_Map& wavefunctionMap)
{
	return Epetra_Vector_Ptr(new Epetra_Vector(View, wavefunctionMap, (double*)psi->GetData().data()));
}
template Epetra_Vector_Ptr CreateWavefunctionEpetraView<1>(Wavefunction<1>::Ptr psi, Epetra_Map& wavefunctionMap);
template Epetra_Vector_Ptr CreateWavefunctionEpetraView<2>(Wavefunction<2>::Ptr psi, Epetra_Map& wavefunctionMap);
template Epetra_Vector_Ptr CreateWavefunctionEpetraView<3>(Wavefunction<3>::Ptr psi, Epetra_Map& wavefunctionMap);
template Epetra_Vector_Ptr CreateWavefunctionEpetraView<4>(Wavefunction<4>::Ptr psi, Epetra_Map& wavefunctionMap);

/*
 * Creates a distributed Epetra_FECrsMatrix from a tensor potential
 */
template<int Rank> 
Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix(typename Wavefunction<Rank>::Ptr psi, blitz::Array<cplx, Rank> potentialData, list pyLocalBasisPairs, double cutoff)
{
	//Create Comm object
	typename DistributedModel<Rank>::Ptr distr = psi->GetRepresentation()->GetDistributedModel();
	Epetra_Comm_Ptr comm = CreateDistributedModelEpetraComm<Rank>(distr);

	//Create Matrix Structures
	int globalSize = blitz::product( psi->GetRepresentation()->GetFullShape() );
	Epetra_Map processorMap(globalSize*2, 0, *comm);

	//Get strides of the global Rank-dimensional array
	blitz::TinyVector<int, Rank> globalShape = psi->GetRepresentation()->GetFullShape();
	blitz::TinyVector<int, Rank> globalStrides;
	globalStrides(Rank-1) = 1;
	for (int rank=Rank-2; rank>=0; rank--)
	{
		globalStrides(rank) = globalStrides(rank+1)	* globalShape(rank+1);
	}

	return CreateTensorPotentialEpetraMatrix(processorMap, potentialData, pyLocalBasisPairs, globalStrides, cutoff);
}
template Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix<1>(Wavefunction<1>::Ptr psi, blitz::Array<cplx, 1> potentialData, list pyLocalBasisPairs, double cutoff);
template Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix<2>(Wavefunction<2>::Ptr psi, blitz::Array<cplx, 2> potentialData, list pyLocalBasisPairs, double cutoff);
template Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix<3>(Wavefunction<3>::Ptr psi, blitz::Array<cplx, 3> potentialData, list pyLocalBasisPairs, double cutoff);
template Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix<4>(Wavefunction<4>::Ptr psi, blitz::Array<cplx, 4> potentialData, list pyLocalBasisPairs, double cutoff);


/*
 * Creates a Epetra_FECrsMatrix from a tensor potential
 */
template<int Rank> 
Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix(Epetra_Map &processorMap, blitz::Array<cplx, Rank> potentialData, list pyLocalBasisPairs, blitz::TinyVector<int, Rank> globalStrides, double cutoff, bool globalAssemble)
{
	Epetra_FECrsMatrix_Ptr matrix( new Epetra_FECrsMatrix(Copy, processorMap, 0) );

	CopyTensorPotentialToEpetraMatrix(matrix, potentialData, pyLocalBasisPairs, globalStrides, cutoff);

	if (globalAssemble)
	{
		matrix->GlobalAssemble();
	}

	return matrix;
}
template Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix(Epetra_Map &processorMap, blitz::Array<cplx, 1> potentialData, list pyLocalBasisPairs, blitz::TinyVector<int, 1> globalStrides, double cutoff, bool globalAssemble);
template Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix(Epetra_Map &processorMap, blitz::Array<cplx, 2> potentialData, list pyLocalBasisPairs, blitz::TinyVector<int, 2> globalStrides, double cutoff, bool globalAssemble);
template Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix(Epetra_Map &processorMap, blitz::Array<cplx, 3> potentialData, list pyLocalBasisPairs, blitz::TinyVector<int, 3> globalStrides, double cutoff, bool globalAssemble);
template Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix(Epetra_Map &processorMap, blitz::Array<cplx, 4> potentialData, list pyLocalBasisPairs, blitz::TinyVector<int, 4> globalStrides, double cutoff, bool globalAssemble);

/*
 * Copy data from tensor potential to already existing Epetra matrix
 */
template<int Rank> 
void CopyTensorPotentialToEpetraMatrix(Epetra_FECrsMatrix_Ptr epetraMatrix, blitz::Array<cplx, Rank> potentialData, list pyLocalBasisPairs, blitz::TinyVector<int, Rank> globalStrides, double cutoff)
{
	blitz::Array<int, 2> indexArray;
	indexArray.resize(4 * potentialData.size(), 2);
	indexArray = -100;

	double sqrCutoff = sqr(cutoff);

	//Setup structures for calculating matrix row/col indices from the 
	//basis pairs in the tensor potential
	blitz::TinyVector< blitz::Array<int, 2>, Rank > localBasisPairs;
	for (int rank=0; rank<Rank; rank++)
	{
		localBasisPairs(rank).reference( boost::python::extract< blitz::Array<int, 2> >(pyLocalBasisPairs[rank]) );
	}

	//Iterate over all items in potentialData
	typename blitz::Array<cplx, Rank>::iterator it = potentialData.begin();
	for (int linearCount=0; linearCount<potentialData.size(); linearCount++)
	{
		int globalRow = 0;
		int globalCol = 0;
		for (int rank=0; rank<Rank; rank++)
		{
			int rankPos = it.position()(rank);
			globalRow += globalStrides(rank) * localBasisPairs(rank)(rankPos, 0);
			globalCol += globalStrides(rank) * localBasisPairs(rank)(rankPos, 1);
		}

		double realVal = real(*it);
		double imagVal = imag(*it);
		
		/*
		 * Because epetra does not support complex natively,
		 * each matrix element is a 2x2 block
		 *
		 * (A_r  -A_i )  (c_r)  =  (A_r + i A_i) * (c_r + i c_i)  =  A * c
		 * (A_i   A_r )  (c_i)
		 *
		 * Detect if A_i or A_r is zero to avoid redundant elements
		 */

		//Insert values into matrix
		if (sqr(realVal) > sqrCutoff)
		{
			int r,c;
			r = 2*globalRow; c = 2*globalCol;
			indexArray(4*linearCount, 0) = r;
			indexArray(4*linearCount, 1) = c;
			epetraMatrix->InsertGlobalValues(r, 1, &realVal, &c);
			r++; c++;
			epetraMatrix->InsertGlobalValues(r, 1, &realVal, &c);
			indexArray(4*linearCount+1, 0) = r;
			indexArray(4*linearCount+1, 1) = c;
		}
		if (sqr(imagVal) > sqrCutoff)
		{
			int r,c;
			r = 2*globalRow; c = 2*globalCol+1;
			indexArray(4*linearCount+2, 0) = r;
			indexArray(4*linearCount+2, 1) = c;
			epetraMatrix->InsertGlobalValues(r, 1, &imagVal, &c);
			imagVal = -imagVal;
			epetraMatrix->InsertGlobalValues(c, 1, &imagVal, &r);
			indexArray(4*linearCount+3, 0) = c;
			indexArray(4*linearCount+3, 1) = r;
		}
		++it;	
	}
}
template void CopyTensorPotentialToEpetraMatrix(Epetra_FECrsMatrix_Ptr epetraMatrix, blitz::Array<cplx, 1> potentialData, list pyLocalBasisPairs, blitz::TinyVector<int, 1> globalStrides, double cutoff);
template void CopyTensorPotentialToEpetraMatrix(Epetra_FECrsMatrix_Ptr epetraMatrix, blitz::Array<cplx, 2> potentialData, list pyLocalBasisPairs, blitz::TinyVector<int, 2> globalStrides, double cutoff);
template void CopyTensorPotentialToEpetraMatrix(Epetra_FECrsMatrix_Ptr epetraMatrix, blitz::Array<cplx, 3> potentialData, list pyLocalBasisPairs, blitz::TinyVector<int, 3> globalStrides, double cutoff);
template void CopyTensorPotentialToEpetraMatrix(Epetra_FECrsMatrix_Ptr epetraMatrix, blitz::Array<cplx, 4> potentialData, list pyLocalBasisPairs, blitz::TinyVector<int, 4> globalStrides, double cutoff);


#endif // use trilinos


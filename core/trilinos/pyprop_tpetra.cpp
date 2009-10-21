#ifdef PYPROP_USE_TRILINOS
#include "pyprop_tpetra.h"

#include <Teuchos_OpaqueWrapper.hpp>

#ifdef EPETRA_MPI
#include <mpi.h>
#include <Teuchos_DefaultMpiComm.hpp>
#else
#include <Teuchos_DefaultSerialComm.hpp>
#endif



/*
 * Creates an Tpetra_Comm object from a DistributedModel. It should
 * use the distributed model to get an MPI_Communicator, but right now
 * DistributedModel uses MPI_COMM_WORLd.
 */
template<int Rank>
Tpetra_Comm_Ptr CreateDistributedModelTpetraComm(typename DistributedModel<Rank>::Ptr distr)
{
	using namespace Teuchos;
#ifdef EPETRA_MPI
	//should get MPI_COMM from distributed model
	return RCP< MpiComm<int> >( new MpiComm<int>(opaqueWrapper(MPI_COMM_WORLD)) );
#else
	return RCP< SerialComm<int> >( new SerialComm<int>() );
#endif
}

template Tpetra_Comm_Ptr CreateDistributedModelTpetraComm<1>(DistributedModel<1>::Ptr distr);
template Tpetra_Comm_Ptr CreateDistributedModelTpetraComm<2>(DistributedModel<2>::Ptr distr);
template Tpetra_Comm_Ptr CreateDistributedModelTpetraComm<3>(DistributedModel<3>::Ptr distr);
template Tpetra_Comm_Ptr CreateDistributedModelTpetraComm<4>(DistributedModel<4>::Ptr distr);



/*
 * Creates an Tpetra_Map corresponding to how a DistributedModel
 * distributes a wavefunction.
 *
 * This map can be used to create Tpetra_Vectors that are array-compatible
 * with Wavefunctions, i.e. has the same memory layout and distribution
 */
template<int Rank>
Tpetra_Map_Ptr CreateWavefunctionTpetraMap(typename Wavefunction<Rank>::Ptr psi)
{
	typedef blitz::Array<cplx, Rank> ArrayType;
	

	ArrayType data = psi->GetData();
	typename DistributedModel<Rank>::Ptr distr = psi->GetRepresentation()->GetDistributedModel();
	Tpetra_Comm_Ptr comm = CreateDistributedModelTpetraComm<Rank>(distr);

	int localSize = data.size();
	int* myElems = new int[localSize];

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

		myElems[linearCount] = globalPos;
		++it;	
	}
	
	Teuchos::ArrayView<int> myElemsView(myElems, localSize);
	return Tpetra_Map_Ptr( new Tpetra_Map(-1, myElemsView, 0, comm) );
}
template Tpetra_Map_Ptr CreateWavefunctionTpetraMap<1>(Wavefunction<1>::Ptr psi);
template Tpetra_Map_Ptr CreateWavefunctionTpetraMap<2>(Wavefunction<2>::Ptr psi);
template Tpetra_Map_Ptr CreateWavefunctionTpetraMap<3>(Wavefunction<3>::Ptr psi);
template Tpetra_Map_Ptr CreateWavefunctionTpetraMap<4>(Wavefunction<4>::Ptr psi);

/*
 * Creates an Tpetra_Vector that is a view into a wavefunction data.
 *   - the wavefunctionMap should be created with CreateWavefunctionTpetraMap(...)
 */
template<int Rank>
Tpetra_Vector_Ptr CreateWavefunctionTpetraView(typename Wavefunction<Rank>::Ptr psi, Tpetra_Map& wavefunctionMap)
{
	Teuchos::ArrayView<cplx> myElemsView((cplx*)psi->GetData().data(), psi->GetData().size());
	return Tpetra_Vector_Ptr(new Tpetra_Vector(wavefunctionMap, myElemsView));
}
template Tpetra_Vector_Ptr CreateWavefunctionTpetraView<1>(Wavefunction<1>::Ptr psi, Tpetra_Map& wavefunctionMap);
template Tpetra_Vector_Ptr CreateWavefunctionTpetraView<2>(Wavefunction<2>::Ptr psi, Tpetra_Map& wavefunctionMap);
template Tpetra_Vector_Ptr CreateWavefunctionTpetraView<3>(Wavefunction<3>::Ptr psi, Tpetra_Map& wavefunctionMap);
template Tpetra_Vector_Ptr CreateWavefunctionTpetraView<4>(Wavefunction<4>::Ptr psi, Tpetra_Map& wavefunctionMap);


#endif // use trilinos


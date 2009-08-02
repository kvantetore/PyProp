#ifndef PYPROP_EPETRA
#define PYPROP_EPETRA
#ifdef PYPROP_USE_TRILINOS

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/representation.h"

#include <Epetra_Comm.h>
#include <Epetra_Time.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_SerialComm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>

#ifndef SINGLEPROC
#include <mpi.h>
#include <Epetra_MpiComm.h>
#endif

/*
 * 
 * Functions to Create Epetra-objects from pyprop objects
 * see pyprop_epetra.cpp for details
 *
 */

typedef shared_ptr<Epetra_Comm> Epetra_Comm_Ptr;
typedef shared_ptr<Epetra_Vector> Epetra_Vector_Ptr;
typedef shared_ptr<Epetra_FECrsMatrix> Epetra_FECrsMatrix_Ptr;
typedef shared_ptr<Epetra_CrsMatrix> Epetra_CrsMatrix_Ptr;
typedef shared_ptr<Epetra_Map> Epetra_Map_Ptr;

template<int Rank>
Epetra_Comm_Ptr CreateDistributedModelEpetraComm(typename DistributedModel<Rank>::Ptr distr);

template<int Rank>
Epetra_Map_Ptr CreateWavefunctionEpetraMap(typename Wavefunction<Rank>::Ptr psi);


template<int Rank>
Epetra_Vector_Ptr CreateWavefunctionEpetraView(typename Wavefunction<Rank>::Ptr psi, Epetra_Map& wavefunctionMap);

template<int Rank> 
Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix(typename Wavefunction<Rank>::Ptr psi, blitz::Array<cplx, Rank> potentialData, boost::python::list pyLocalBasisPairs, double cutoff);

template<int Rank> 
Epetra_FECrsMatrix_Ptr CreateTensorPotentialEpetraMatrix(Epetra_Map &processorMap, blitz::Array<cplx, Rank> potentialData, boost::python::list pyLocalBasisPairs, blitz::TinyVector<int, Rank> globalStrides, double cutoff);


#endif // use trilinos
#endif

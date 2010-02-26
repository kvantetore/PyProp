#ifndef PYPROP_TPETRA
#define PYPROP_TPETRA
#ifdef PYPROP_USE_TRILINOS
#ifdef PYPROP_USE_TRILINOS_TPETRA

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/representation.h"

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Operator.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_Comm.hpp>

/*
 * 
 * Functions to Create Epetra-objects from pyprop objects
 * see pyprop_epetra.cpp for details
 *
 */

typedef Tpetra::Vector<cplx> Tpetra_Vector;
typedef Tpetra::MultiVector<cplx> Tpetra_MultiVector;
typedef Tpetra::Map<int> Tpetra_Map;
typedef Teuchos::Comm<int> Tpetra_Comm;
typedef Tpetra::Operator<cplx> Tpetra_Operator;

typedef Teuchos::RCP<Tpetra_Comm> Tpetra_Comm_Ptr;
typedef shared_ptr<Tpetra_Vector> Tpetra_Vector_Ptr;
typedef shared_ptr<Tpetra_Map> Tpetra_Map_Ptr;

template<int Rank>
Tpetra_Comm_Ptr CreateDistributedModelTpetraComm(typename DistributedModel<Rank>::Ptr distr);

template<int Rank>
Tpetra_Map_Ptr CreateWavefunctionTpetraMap(typename Wavefunction<Rank>::Ptr psi);

template<int Rank>
Tpetra_Vector_Ptr CreateWavefunctionTpetraView(typename Wavefunction<Rank>::Ptr psi, Tpetra_Map& wavefunctionMap);

#endif // use tpetra
#endif // use trilinos
#endif

#ifndef DISTRIBUTEDOVERLAPMATRIX_H
#define DISTRIBUTEDOVERLAPMATRIX_H

#include "../wavefunction.h"
#include <boost/shared_ptr.hpp>
#include <Epetra_MultiVector.h>
#include "../trilinos/pyprop_epetra.h"

#include "distributedoverlapmatrix.h"

template<int Rank>
class DistributedOverlapMatrix
{
public:
	typedef shared_ptr< DistributedOverlapMatrix<Rank> > Ptr;

private:
	typedef shared_ptr<Epetra_MultiVector> Epetra_MultiVector_Ptr;
	typedef blitz::TinyVector<Epetra_CrsMatrix_Ptr, Rank> CrsVector;
	typedef blitz::TinyVector<Epetra_Map_Ptr, Rank> EpetraMapVector;
	
	CrsVector OverlapMatrices; 
	EpetraMapVector WavefunctionMaps;
	blitz::TinyVector<bool, Rank> IsSetupRank;
	blitz::TinyVector<Epetra_MultiVector_Ptr, Rank> InputVector;
	blitz::TinyVector<Epetra_MultiVector_Ptr, Rank> OutputVector;

	virtual void SetupOverlapMatrixRank(Wavefunction<Rank> &srcPsi, int opRank);

public:

	DistributedOverlapMatrix() 
	{
		IsSetupRank = false;
	}
	virtual ~DistributedOverlapMatrix() {}
	
	virtual void SetupRank(Wavefunction<Rank> &srcPsi, int opRank);
	virtual Epetra_MultiVector_Ptr SetupMultivector(Wavefunction<Rank> &srcPsi, int opRank);
	virtual void MultiVectorToWavefunction(Wavefunction<Rank> &psi, Epetra_MultiVector_Ptr vec, int opRank);
	virtual void WavefunctionToMultiVector(Wavefunction<Rank> &psi, Epetra_MultiVector_Ptr vec, int opRank);
	virtual void MultiplyOverlapRank(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &destPsi, int opRank, bool fastAlgo);
	virtual void MultiplyOverlapRank(Wavefunction<Rank> &psi, int opRank, bool fastAlgo)
	{
		MultiplyOverlapRank(psi, psi, opRank, fastAlgo);
	}
	virtual void SolveOverlapRank(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &destPsi, int opRank, bool fastAlgo);
	virtual void SolveOverlapRank(Wavefunction<Rank> &psi, int opRank, bool fastAlgo)
	{
		SolveOverlapRank(psi, psi, opRank, fastAlgo);
	}
};

#endif


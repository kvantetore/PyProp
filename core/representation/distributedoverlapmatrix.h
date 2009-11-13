#ifndef DISTRIBUTEDOVERLAPMATRIX_H
#define DISTRIBUTEDOVERLAPMATRIX_H

#include "../wavefunction.h"
#include <boost/shared_ptr.hpp>
#include <Epetra_MultiVector.h>
#include "../trilinos/pyprop_epetra.h"

template<int Rank>
class DistributedOverlapMatrix
{

private:
	typedef shared_ptr<Epetra_MultiVector> Epetra_MultiVector_Ptr;
	typedef blitz::TinyVector<Epetra_CrsMatrix_Ptr, Rank> CrsVector;
	typedef blitz::TinyVector<Epetra_Map_Ptr, Rank> EpetraMapVector;
	
	CrsVector OverlapMatrices; 
	EpetraMapVector WavefunctionMaps;

public:
	DistributedOverlapMatrix() {}
	virtual ~DistributedOverlapMatrix() {}
	
	virtual void SetupRank(Wavefunction<Rank> &srcPsi, int opRank);
	virtual Epetra_MultiVector_Ptr SetupMultivector(Wavefunction<Rank> &srcPsi, int opRank);
	virtual void MultiVectorToWavefunction(Wavefunction<Rank> &psi, Epetra_MultiVector_Ptr vec, int opRank);
	virtual void MultiplyOverlapRank(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &destPsi, int opRank);
	virtual void SolveOverlapRank(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &destPsi, int opRank);
};

#endif


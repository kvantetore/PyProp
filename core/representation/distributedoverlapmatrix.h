#ifndef DISTRIBUTEDOVERLAPMATRIX_H
#define DISTRIBUTEDOVERLAPMATRIX_H

#include "../wavefunction.h"
#include <boost/shared_ptr.hpp>
#include <Epetra_MultiVector.h>
#include "../trilinos/pyprop_epetra.h"

#include "distributedoverlapmatrix.h"

#include <Epetra_LinearProblem.h>

#include <Amesos_ConfigDefs.h>
#include <Amesos.h>

template<int Rank>
class DistributedOverlapMatrix
{
public:
	typedef shared_ptr< DistributedOverlapMatrix<Rank> > Ptr;

private:
	typedef shared_ptr<Epetra_MultiVector> Epetra_MultiVector_Ptr;
	typedef blitz::TinyVector<Epetra_CrsMatrix_Ptr, Rank> CrsVector;
	typedef blitz::TinyVector<Epetra_Map_Ptr, Rank> EpetraMapVector;
	typedef shared_ptr<Epetra_LinearProblem> Epetra_LinearProblem_Ptr;
	typedef shared_ptr<Amesos> Amesos_Ptr;
	typedef shared_ptr<Amesos_BaseSolver> Amesos_BaseSolver_Ptr;

	bool HasPsi;
	CrsVector OverlapMatrices; 
	EpetraMapVector WavefunctionMaps;
	blitz::TinyVector<bool, Rank> IsSetupRank;
	blitz::TinyVector<Epetra_MultiVector_Ptr, Rank> InputVector;
	blitz::TinyVector<Epetra_MultiVector_Ptr, Rank> OutputVector;

	typename Wavefunction<Rank>::Ptr Psi;

	blitz::Array<cplx, Rank> TempData;
	blitz::Array<cplx, Rank> InputData;
	blitz::Array<cplx, Rank> OutputData;

	//Amesos objects for each rank
	blitz::TinyVector<Epetra_LinearProblem_Ptr, Rank> EpetraProblems;
	blitz::TinyVector<Amesos_BaseSolver_Ptr, Rank> Solvers;
	blitz::TinyVector<Amesos_Ptr, Rank> Factories;

	virtual void SetupOverlapMatrixRank(Wavefunction<Rank> &srcPsi, int opRank);
	virtual void SetupOverlapSolvers(Wavefunction<Rank> &srcPsi, int opRank);

public:

	DistributedOverlapMatrix() 
	{
		IsSetupRank = false;
		HasPsi = false;
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


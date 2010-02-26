#ifndef DISTRIBUTEDMODEL_H
#define DISTRIBUTEDMODEL_H

#include "../common.h"
#include "../wavefunction.h"
#include "distribution.h"

#include <mpi.h>

//Forward declaration of arraytranspose
template<int Rank> class ArrayTranspose;

template<int Rank>
class DistributedModel
{
public:
	typedef shared_ptr< DistributedModel<Rank> > Ptr;
	typedef ArrayTranspose<Rank> TransposeType;
	typedef shared_ptr<TransposeType> TransposePtr;
	typedef shared_ptr< DistributedModel<1> > DistributedModel1DPtr;

private:
	Distribution::Ptr CurrentDistribution;
	TransposePtr Transpose;

	static bool MPIDisabled;
	
public:

	int ProcId;
	int ProcCount;
	
	//constructor
	DistributedModel();
	DistributedModel(DistributionPtr distrib);
	DistributedModel(const DistributedModel &other)
	{
		this->CurrentDistribution = DistributionPtr(new Distribution(*other.CurrentDistribution));
		this->Transpose = other.Transpose;
		this->ProcCount = other.ProcCount;
		this->ProcId = other.ProcId;
	}

	/*
	 * Create a 1D distributed model sharing CurrentDistribution with this DistributedModel
	 * This is used by CombinedRepresentation to construct lower dimensional representations which
	 * still knows the correct distribution
	 */
	DistributedModel1DPtr CreateSubDistributedModel()
	{
		return DistributedModel1DPtr( new DistributedModel<1>(CurrentDistribution) );
	}

	bool IsSingleProc()
	{
		return ProcCount == 1;
	}

	bool IsFirstProc()
	{
		return ProcId == 0;
	}

	bool IsLastProc()
	{
		return ProcId == ProcCount -1;
	}
	
	Distribution::DataArray GetDistribution()
	{
		return CurrentDistribution->GetDistribution();
	}

	void SetDistribution(Distribution::DataArray newDistrib)
	{
		return CurrentDistribution->SetDistribution(newDistrib);
	}


/*
	template<int Rank2>
	void ShareDistribution(DistributedModel<Rank2>::Ptr other)
	{
		this->CurrentDistribution = other->CurrentDistribution;
	}
*/

	TransposePtr GetTranspose()
	{
		return Transpose;
	}
	
	template<class T> const blitz::Array<T, 1> GetLocalArray(const blitz::Array<T, 1> &array, int rank)
	{
		return array(GetLocalIndexRange(array.extent(0), rank));
	}

	static void ForceSingleProc()
	{
		MPIDisabled = true;
	}

	void SetupMPI();
	blitz::TinyVector<int, Rank> CreateInitialShape(const blitz::TinyVector<int, Rank> &fullShape);
	blitz::TinyVector<int, Rank> GetGlobalShape(const blitz::TinyVector<int, Rank> localShape);
	int GetLocalStartIndex(int globalSize, int currentRank);
	blitz::Range GetLocalIndexRange(int globalSize, int currentRank);
	void ApplyConfigSection(const ConfigSection &cfg);

	int GetLocalStartIndex(int globalSize, int currentRank, int procId);
	blitz::Range GetLocalIndexRange(int globalSize, int currentRank, int procId);
	

	double GetGlobalSum(double localValue);
	cplx GetGlobalSum(cplx localValue);
	void GetGlobalSum(blitz::Array<cplx, 1> &in, blitz::Array<cplx, 1> &out);
	void GlobalBarrier();

	bool IsDistributedRank(int rank);
	void ChangeDistribution(Wavefunction<Rank> &psi, const Distribution::DataArray &newDistrib, int destBufferName);

	MPI_Comm GetGroupCommRank(int rank);

	static void InitMPI(int argc, char* argv[]);
	static void FinalizeMPI();
};

template<int Rank>
bool DistributedModel<Rank>::MPIDisabled = false;

#endif


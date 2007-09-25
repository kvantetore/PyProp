
#include "distributedmodel.h"
#include "blitztranspose.h"
#include "../representation/representation.h"

#include <mpi.h>

using namespace blitz;

template<int Rank>
void DistributedModel<Rank>::InitMPI(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
}

template<int Rank>
void DistributedModel<Rank>::FinalizeMPI()
{
	MPI_Finalize();
}

template<int Rank>
DistributedModel<Rank>::DistributedModel() 
{
	SetupMPI();
}

template<int Rank>
DistributedModel<Rank>::DistributedModel(DistributionPtr distrib)
{
	SetupMPI();

	int procRank = distrib->GetProcRank();
	if (!IsSingleProc())
	{
		Transpose = TransposePtr( new ArrayTranspose<Rank>(procRank));
	}
	CurrentDistribution = distrib;
}

template<int Rank>
void DistributedModel<Rank>::SetupMPI()
{
	if (MPIDisabled)
	{
		this->ProcId = 0;
		this->ProcCount = 1;
		std::cout << "MPI Disabled" << endl;
	}
	else
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &this->ProcId);
		MPI_Comm_size(MPI_COMM_WORLD, &this->ProcCount);
		std::cout << "Proc " << (1+ProcId) << "/" << ProcCount << std::endl;
	}
}


template<int Rank>
void DistributedModel<Rank>::ApplyConfigSection(const ConfigSection &cfg)
{
	//How many dimensions in the processor array
	int procRank = -1;
	cfg.Get("proc_array_rank", procRank);
	if (!IsSingleProc())
	{
		Transpose = TransposePtr( new TransposeType(procRank) );
	}

	//How is the initial distribution
	CurrentDistribution = DistributionPtr( new Distribution(procRank) );
	Distribution::DataArray distr(procRank);
	cfg.Get("initial_distribution", distr);
	CurrentDistribution->SetDistribution(distr);
}

/**

*/
template<int Rank>
blitz::TinyVector<int, Rank> DistributedModel<Rank>::CreateInitialShape(const blitz::TinyVector<int, Rank> &fullShape)
{
	if (IsSingleProc())
	{
		return fullShape;
	}

	return Transpose->CreateDistributedShape(fullShape, CurrentDistribution->GetDistribution());
}


template<int Rank>
bool DistributedModel<Rank>::IsDistributedRank(int rank)
{
	if (IsSingleProc())
	{
		return false;
	}

	Distribution::DataArray distrib(CurrentDistribution->GetDistribution());
	for (int i=0; i<CurrentDistribution->GetProcRank(); i++)
	{
		if (rank == distrib(i))
		{
			int Np = Transpose->GetProcGridShape()(i);
			if (Np > 1)
			{
				return true;
			}
		}
	}
	return false;
}


/* Performs the necessary communication with the other processes to change the distribution
 * of the wavefunction
 *    - psi             The wavefunction
 *    - newDistrib:     The new distribution. This should differ from current distribution only by
 *                      one rank currently distributed
 *    - destBufferName: The buffer in psi which should be used to recieve the new distribution.
 */
template<int Rank>
void DistributedModel<Rank>::ChangeDistribution(Wavefunction<Rank> &psi, const Distribution::DataArray &newDistrib, int destBufferName)
{
	if (IsSingleProc())
	{
		throw std::runtime_error("Error! Should not call ChangeDistribution if we are singleproc.");
	}

	typename Wavefunction<Rank>::DataArray src(psi.GetData());
	typename Wavefunction<Rank>::DataArray dst(psi.GetData(destBufferName));
	typename Wavefunction<Rank>::IndexVector fullShape(psi.GetRepresentation().GetFullShape());

	Transpose->Transpose(fullShape, src, CurrentDistribution->GetDistribution(), dst, newDistrib);

	psi.SetActiveBuffer(destBufferName);
	CurrentDistribution->SetDistribution(newDistrib);
}

template<>
void DistributedModel<1>::ChangeDistribution(Wavefunction<1> &psi, const Distribution::DataArray &newDistrib, int destBufferName)
{
	cout << "Error! This is wrong. ChangeDistribution should never be called on a 1D DistributedModel. "
	     << "Most likely, you are calling ChangeDistribution on a SubDistribution. " << endl;
}

/**
 * Returns the index in the virtual super wavefunction
 * on which this proc starts (on the rank currentRank)
 */
template<int Rank>
int DistributedModel<Rank>::GetLocalStartIndex(int globalSize, int currentRank)
{
	if (IsSingleProc())
	{
		return 0;
	}

	Distribution::DataArray distrib(CurrentDistribution->GetDistribution());
	for (int i=0; i<CurrentDistribution->GetProcRank(); i++)
	{
		if (distrib(i) == currentRank)
		{
			return Transpose->GetLocalStartIndex(globalSize, i);
		}
	}
	return 0;
}

/** 
 * Returns the local portion of the range of currentRank in 
 * the global virtual super-wavefunction.  
 */	
template<int Rank>
blitz::Range DistributedModel<Rank>::GetLocalIndexRange(int globalSize, int currentRank)
{
	if (IsSingleProc())
	{
		return blitz::Range(0, globalSize-1);
	}

	Distribution::DataArray distrib(CurrentDistribution->GetDistribution());
	for (int i=0; i<CurrentDistribution->GetProcRank(); i++)
	{
		if (distrib(i) == currentRank)
		{
			return Transpose->GetLocalRange(globalSize, i);
		}
	}
	return blitz::Range(0, globalSize-1);
}


template<int Rank>
double DistributedModel<Rank>::GetGlobalSum(double localValue)
{
	if (IsSingleProc())
	{
		return localValue;
	}

	double globalValue = 0;
	MPI_Allreduce(&localValue, &globalValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalValue;
}

template<int Rank>
cplx DistributedModel<Rank>::GetGlobalSum(cplx localValue)
{
	if (IsSingleProc())
	{
		return localValue;
	}

	cplx globalValue = 0;
	MPI_Allreduce(&localValue, &globalValue, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalValue;
}

template<int Rank>
void DistributedModel<Rank>::GlobalBarrier()
{
	if (!IsSingleProc())
	{
		MPI_Barrier(MPI_COMM_WORLD);
	}
}



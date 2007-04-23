
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
	Transpose = TransposePtr( new ArrayTranspose<Rank>(procRank));
	CurrentDistribution = distrib;
}

template<int Rank>
void DistributedModel<Rank>::SetupMPI()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &this->ProcId);
	MPI_Comm_size(MPI_COMM_WORLD, &this->ProcCount);
	std::cout << "Proc " << (1+ProcId) << "/" << ProcCount << std::endl;
}


template<int Rank>
void DistributedModel<Rank>::ApplyConfigSection(const ConfigSection &cfg)
{
	int procRank = 1;
	if (cfg.HasValue("proc_array_rank"))
	{
		cfg.Get("proc_array_rank", procRank);
		cout << "Missing proc_array_rank from Distribution, using default." << endl;
	}
	Transpose = TransposePtr( new TransposeType(procRank) );
	
	CurrentDistribution = DistributionPtr( new Distribution(procRank) );
	Distribution::DataArray distr(procRank);
	if (cfg.HasValue("initial_distribution"))
	{
		cfg.Get("initial_distribution", distr);
	}
	else
	{
		//default is to initially distribute the procRank first ranks.
		distr = blitz::tensor::i;
		cout << "Missing initial_distribution from Distribution, using default." << endl;
	}
	CurrentDistribution->SetDistribution(distr);
}

/**

*/
template<int Rank>
blitz::TinyVector<int, Rank> DistributedModel<Rank>::CreateInitialShape(const blitz::TinyVector<int, Rank> &fullShape)
{
	cout << "Creating distributed shape to " << fullShape;
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
	double globalValue = 0;
	MPI_Allreduce(&localValue, &globalValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalValue;
}



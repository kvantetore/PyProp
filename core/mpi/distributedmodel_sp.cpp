
#include "distributedmodel.h"
#include "../representation/representation.h"

using namespace blitz;

template<int Rank>
void DistributedModel<Rank>::InitMPI(int argc, char* argv[])
{
	cout << "MPI Not enabled, recompile pyprop without -DSINGLEPROC in order to enale MPI" << endl;
}

template<int Rank>
void DistributedModel<Rank>::FinalizeMPI()
{
	cout << "MPI Not enabled, recompile pyprop without -DSINGLEPROC in order to enale MPI" << endl;
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

	CurrentDistribution = distrib;
}

template<int Rank>
void DistributedModel<Rank>::SetupMPI()
{
	this->ProcId = 0;
	this->ProcCount = 1;
}


template<int Rank>
void DistributedModel<Rank>::ApplyConfigSection(const ConfigSection &cfg)
{
}

/**

*/
template<int Rank>
blitz::TinyVector<int, Rank> DistributedModel<Rank>::CreateInitialShape(const blitz::TinyVector<int, Rank> &fullShape)
{
	cout << "Creating distributed shape to " << fullShape;
	return fullShape;
}


template<int Rank>
bool DistributedModel<Rank>::IsDistributedRank(int rank)
{
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
	return 0;
}

/** 
 * Returns the local portion of the range of currentRank in 
 * the global virtual super-wavefunction.  
 */	
template<int Rank>
blitz::Range DistributedModel<Rank>::GetLocalIndexRange(int globalSize, int currentRank)
{
	return blitz::Range(0, globalSize-1);
}


template<int Rank>
double DistributedModel<Rank>::GetGlobalSum(double localValue)
{
	return localValue;
}



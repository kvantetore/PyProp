#ifndef DISTRIBUTEDMODEL_H
#define DISTRIBUTEDMODEL_H

#include "../wavefunction.h"
#include "../configuration.h"
#include "../hltimer.h"
#include <iostream>

enum TransposeDirection
{
	TRANSPOSE_FORWARD = -1,
	TRANSPOSE_BACKWARD = 1
};

enum TransposeModel
{
	TRANSPOSE_CYCLIC 	= 0, 
	TRANSPOSE_ALTERNATE	= 1,
	TRANSPOSE_SEMI		= 2	
};

template<int Rank>
class DistributedModel
{
private:
	unsigned int CurrentStep;
	unsigned int DistributedRank;
	
public:
	TransposeModel Model;
	int ProcId;
	int ProcCount;
	
	double MpiTime; 
	double TotalTime;
	HL::Timer TotalTimer;
	HL::Timer MpiTimer;
	
	//constructor
	DistributedModel();
	
	/**
	Change the distributed representation according to Model
	*/
	void ChangeRepresentation(Wavefunction<Rank> &psi);
	
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
	
	int GetInitialDistributedRank()
	{
		//For now, the first rank is distributed initially
		return 0;
	}
	
	bool IsDistributedRank(int rank)
	{
		if (IsSingleProc())
		{
			return false;
		}
		return GetDistributedRank() == rank;
	}
	
	/**
	
	*/
	blitz::TinyVector<int, Rank> CreateInitialShape(const blitz::TinyVector<int, Rank> &fullShape)
	{
		blitz::TinyVector<int, Rank> shape = fullShape;
		int distribRank = GetInitialDistributedRank();
		shape(distribRank) = shape(distribRank) / this->ProcCount;
		return shape;
	}
	
	/**
	Returns the rank which is currently distributed among
	the procs.
	*/
	int GetDistributedRank()
	{
		return DistributedRank;
	}
	
	/**
	Returns if the currently distributed rank has maximum stride.
	*/
	bool HasDistributedRangeMaxStride()
	{
		return GetDistributedRank() == Rank - 1;
	}
	
	/**
	Returns the index for the global virtual super-wavefunction.
	*/
	blitz::Range GetFullIndexRange(const Wavefunction<Rank> &psi, int currentRank)
	{
		int extent = psi.Data.extent(currentRank);
		if (currentRank != GetDistributedRank()) 
		{
			return blitz::Range(0, extent - 1);
		} 
		else
		{
			return blitz::Range(0, extent * ProcCount - 1);
		}
	}
	
	/**
	Returns the index of the local part of the wavefunction
	in current-proc coordinates.
	*/
	blitz::Range GetLocalIndexRange(const Wavefunction<Rank> &psi, int currentRank)
	{
		int extent = psi.Data.extent(currentRank);
		return blitz::Range(0, extent - 1);
	}
	

	/**
	 * Returns the index in the virtual super wavefunction
	 * on which this proc starts
	 */
	int GetGlobalStartIndex(int globalSize, int currentRank)
	{
		int base = 0;
		if (currentRank == GetDistributedRank()) 
		{
			base = globalSize * ProcId;
		}
		return base;
	}
	
	/** 
	 * Returns the local portion of the range of currentRank in 
	 * the global virtual super-wavefunction. i.e. for a local
	 * rank it returns the same as GetLocalIndex, and for a 
	 * distributed rank, it returns (extent*ProcId -> extent * (procid + 1))
	 */	
	blitz::Range GetGlobalIndexRange(int globalSize, int currentRank)
	{
		int base = GetGlobalStartIndex(globalSize, currentRank);
		return blitz::Range(base, base + globalSize - 1);
	}

	template<class T> const blitz::Array<T, 1> GetLocalArray(const blitz::Array<T, 1> &array, int rank)
	{
		//Should be implemented correctly, need to set up a proper way to figure out the
		//partitioning, etc.
		if (IsDistributedRank(rank))
		{
			cout << "Warning: Should implement DistributedModel::GetLocalArray properly" << endl;
		}
		return array;
	}


	void ApplyConfigSection(const ConfigSection &cfg);

		
private:
	//temp data	
	blitz::Array<cplx, 2> tempA;
	blitz::Array<cplx, 2> tempB;

	int BufferName1; //Should perhaps be generalized to more than two buffers
	int BufferName2; //but right now, mpi support is pretty broken anyway

	void Transpose(blitz::Array<cplx,Rank> &A, blitz::Array<cplx,Rank> &B, TransposeDirection direction);
	int TransposeSemi(blitz::Array<cplx,Rank> &A, blitz::Array<cplx,Rank> &B, int rank, TransposeDirection direction);
	void Transpose2d(blitz::Array<cplx,2> &A, blitz::Array<cplx,2> &B);
	void Transpose2dBack(blitz::Array<cplx,2> &A, blitz::Array<cplx,2> &B);
	
};



#endif


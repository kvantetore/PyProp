#ifndef REPRESENTATION_H
#define REPRESENTATION_H

#include "../common.h"
#include "../mpi/distributedmodel.h"

template<int Rank> class Wavefunction;

template<int Rank>
class Representation
{
public:
	typedef boost::shared_ptr< DistributedModel<Rank> > DistributedModelPtr;
	
private:
	DistributedModelPtr Distrib;	//Knows how the wavefunction is distributed
	int BaseRank;					//Used when this representation is a part of a combined representation
									//All rank parameters will use the global rank. It is the responsibility
									//of each representation class to translate down to the effective rank
									//(subtract BaseRank)

public:
	//Constructors
	Representation () : BaseRank(0) {}
	virtual ~Representation() {}
	
	inline void SetDistributedModel(DistributedModelPtr distrib)
	{
		Distrib = distrib;
	}
	
	inline DistributedModel<Rank>& GetDistributedModel()
	{
		return *Distrib;
	}
	
	blitz::TinyVector<int, Rank> GetInitialShape() 
	{
		blitz::TinyVector<int, Rank> fullShape = GetFullShape();
		return GetDistributedModel().CreateInitialShape(fullShape);
	}

	int GetBaseRank()
	{
		return BaseRank;
	}

	void SetBaseRank(int baseRank)
	{
		BaseRank = baseRank;
	}

	long int GetId()
	{
		return (long int)this;
	}

	virtual blitz::Array<double, 1> GetLocalGrid(int rank)
	{
		return this->GetDistributedModel().GetLocalArray(GetGlobalGrid(rank), rank);
	}

	//Must override
	virtual blitz::TinyVector<int, Rank> GetFullShape() = 0;
	virtual cplx InnerProduct(const Wavefunction<Rank> &w1, const Wavefunction<Rank> &w2) = 0;
	virtual blitz::Array<double, 1> GetLocalWeights(int rank) = 0;
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank) = 0;
	virtual void ApplyConfigSection(const ConfigSection &config) = 0;

};

#endif

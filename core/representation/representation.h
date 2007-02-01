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
	DistributedModelPtr Distrib;

public:
	//Constructors
	Representation () {}
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

	long int GetId()
	{
		return (long int)this;
	}

	//Must override
	virtual blitz::TinyVector<int, Rank> GetFullShape() = 0;
	virtual cplx InnerProduct(const Wavefunction<Rank> &w1, const Wavefunction<Rank> &w2) = 0;
	virtual blitz::Array<double, 1> GetLocalGrid(const Wavefunction<Rank> &psi, int rank) = 0;
	virtual void ApplyConfigSection(const ConfigSection &config) = 0;

};

#endif

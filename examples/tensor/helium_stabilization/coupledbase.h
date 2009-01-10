#ifndef COUPLEDBASE_H
#define COUPLEDBASE_H

#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>
#include <core/representation/combinedrepresentation.h>
#include <core/representation/coupledspherical/coupledsphericalharmonicrepresentation.h>

template<int Rank>
class CustomPotentialCoupledSphericalBase
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

protected:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank1;
	int RadialRank2;

public:
	CustomPotentialCoupledSphericalBase() {}
	virtual ~CustomPotentialCoupledSphericalBase() {}

	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank1", RadialRank1);
		config.Get("radial_rank2", RadialRank2);
		config.Get("angular_rank", AngularRank);
	}

	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
	{
		if (rank != AngularRank)
		{
			throw std::runtime_error("Only angular rank supports basis pairs");
		}
		AngularBasisPairs.reference(basisPairs.copy());
	}

	virtual BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime) = 0;
};


#endif


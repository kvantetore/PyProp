#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>
#include <core/representation/combinedrepresentation.h>
#include <core/representation/coupledspherical/coupledsphericalharmonicrepresentation.h>

class CoupledSphericalSelectionRule
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

	CoupledSphericalSelectionRule() {}
	virtual ~CoupledSphericalSelectionRule() {}

	virtual bool SelectionRule(CoupledSpherical::CoupledIndex const& left, CoupledSpherical::CoupledIndex const& right) = 0;

	BasisPairList GetBasisPairs(const CoupledSpherical::CoupledSphericalHarmonicRepresentation::Ptr angularRepresentation)
	{
		using namespace CoupledSpherical;

		//Setup angular index pair list
		typedef blitz::TinyVector<int, 2> BasisPair;
		std::vector<BasisPair> list;
		int angCount = angularRepresentation->Range.Count();
		for (int leftIndex=0; leftIndex<angCount; leftIndex++)
		{
			CoupledIndex left = angularRepresentation->Range.GetCoupledIndex(leftIndex);
			
			for (int rightIndex=0; rightIndex<angCount; rightIndex++)
			{
				CoupledIndex right = angularRepresentation->Range.GetCoupledIndex(rightIndex);
				
				if (SelectionRule(left, right))
				{
					list.push_back( BasisPair(leftIndex, rightIndex) );
				}
			}
		}
		//copy list into an array
		int pairIndex = 0;
		blitz::Array<int, 2> indexPairs(list.size(), 2);
		BOOST_FOREACH(BasisPair p, list)
		{
			indexPairs(pairIndex, 0) = p(0);
			indexPairs(pairIndex, 1) = p(1);
			pairIndex++;
		}

		return indexPairs;
	}

};

class CoupledSphericalSelectionRuleR12  : public CoupledSphericalSelectionRule
{
public:
	CoupledSphericalSelectionRuleR12() {}
	virtual ~CoupledSphericalSelectionRuleR12() {}

	virtual bool SelectionRule(CoupledSpherical::CoupledIndex const& left, CoupledSpherical::CoupledIndex const& right)
	{
		return left.L == right.L && left.M == right.M;
	}
};

template<int Rank>
class CustomPotentialEvaluationR12
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;

public:
	CustomPotentialEvaluationR12() {}
	~CustomPotentialEvaluationR12() {}

	void ApplyConfigSection(const ConfigSection &config)
	{
	}

	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
	{
		if (rank != 2)
		{
			throw std::runtime_error("Only rank 2 supports basis pairs");
		}
		AngularBasisPairs.reference(basisPairs.copy());
	}

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == 2)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;

		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CoupledSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< CoupledSphericalHarmonicRepresentation >(repr->GetRepresentation(2));
	
		int r1Count = data.extent(0);
		int r2Count = data.extent(1);
		int angCount = data.extent(2);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(0);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(1);
	
		ClebschGordan cg;

		BasisPairList angBasisPairs = GetBasisPairList(2);

		if (psi->GetRepresentation()->GetDistributedModel()->IsDistributedRank(2)) throw std::runtime_error("Angular rank can not be distributed");
		if (data.extent(0) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(1) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(2) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		data = 0;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
			CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);
	
			if (left.L != right.L) continue;
			if (left.M != right.M) continue;
	
			int L = left.L;
			int M = left.M;
			int l1p = left.l1;
			int l2p = left.l2;
			int l1 = right.l1;
			int l2 = right.l2;

			int minL3 = std::max(std::abs(l1 - l1p), std::abs(l2 - l2p));
			int maxL3 = std::min(l1+l1p, l2+l2p);
	
			for (int l3=minL3; l3<=maxL3; l3++)
			{
				double l3Coeff = 1.0;
				l3Coeff *= CoefficientR12(l1p, l1, l3);
				l3Coeff *= CoefficientR12(l2p, l2, l3);
				l3Coeff *= 1. / (2*l3 + 1.);
				l3Coeff *= cg(l1p, l1, 0, 0, l3, 0);
				l3Coeff *= cg(l2p, l2, 0, 0, l3, 0);

				double l3Sum = 0;
				for (int m1p=-l1; m1p<=l1; m1p++)
				{
					int m2p = M - m1p;
					for (int m1=-right.l1; m1<=l1; m1++)
					{
						int m2 = M - m1;
						for (int m3=-l3; m3<=l3; m3++)
						{
							double cur = 1.0;
							cur *= cg(l1p, l2p, m1p, m2p, L, M);
							cur *= cg(l1, l2, m1, m2, L, M);
							cur *= cg(l1p, l1, -m1p, m1, l3, m3);
							cur *= cg(l2p, l2, -m2p, m2, l3, -m3);
							cur *= std::pow(-1., m1p + m2p + m3);
							l3Sum += cur;
						}
					}
				}
				l3Sum *= l3Coeff; 

				for (int ri1=0; ri1<r1Count; ri1++)
				{
					double r1 = localr1(ri1);
			
					for (int ri2=0; ri2<r2Count; ri2++)
					{
						double r2 = localr2(ri2);
			
						double rmin = std::min(r1, r2);
						double rmax = std::max(r1, r2);
						double rfrac = rmin / rmax;
				
						data(ri1, ri2, angIndex) += l3Sum * std::pow(rfrac, l3) / rmax;
					}
				}
			}
		}
	}
	
	static double CoefficientR12(int l1, int l2, int l3)
	{
		return 	std::sqrt( (2.*l1 + 1.) * (2.*l2 + 1.) / (2.*l3 + 1.) );
	}
};


template<int Rank>
class CoupledSphericalKineticEnergyEvaluator
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;

public:
	CoupledSphericalKineticEnergyEvaluator() {}
	~CoupledSphericalKineticEnergyEvaluator() {}

	double Mass;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", Mass);
	}

	virtual void SetBasisPairs(int rank, const BasisPairList &basisPairs)
	{
		if (rank != 2)
		{
			throw std::runtime_error("Only rank 2 supports basis pairs");
		}
		AngularBasisPairs.reference(basisPairs.copy());
	}

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == 2)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		cout << "Got data of shape " << data.shape() << endl;

		typedef CombinedRepresentation<Rank> CmbRepr;

		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CoupledSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< CoupledSphericalHarmonicRepresentation >(repr->GetRepresentation(2));
	
		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(0);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(1);

		int r1Count = localr1.extent(0);
		int r2Count = localr2.extent(0);
		int angCount = data.extent(2);

		ClebschGordan cg;

		BasisPairList angBasisPairs = GetBasisPairList(2);

		if (psi->GetRepresentation()->GetDistributedModel()->IsDistributedRank(2)) throw std::runtime_error("Angular rank can not be distributed");
		if (data.extent(0) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(1) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(2) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");
	
		for (int ri1=0; ri1<r1Count; ri1++)
		{
			double r1 = localr1(ri1);
	
			for (int ri2=0; ri2<r2Count; ri2++)
			{
				double r2 = localr2(ri2);
	
				for (int angIndex=0; angIndex<angCount; angIndex++)
				{
					int leftIndex = angBasisPairs(angIndex, 0);
					int rightIndex = angBasisPairs(angIndex, 1);
	
					CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
					CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);
	
					if (left.L != right.L) continue;
					if (left.M != right.M) continue;
					if (left.l1 != right.l1) continue;
					if (left.l2 != right.l2) continue;

					double V1 = left.l1 * (left.l1 + 1.) / (2. * Mass * r1 * r1);
					double V2 = left.l2 * (left.l2 + 1.) / (2. * Mass * r2 * r2);
	
					data(ri1, ri2, angIndex) = V1 + V2;
				}
			}
		}
	}
};


template<int Rank>
class CoupledSphericalCoulombPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	double Z;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("z", Z);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r1 = pos(0);
		double r2 = pos(1);
		return - Z/r1 - Z/r2;
	}
};



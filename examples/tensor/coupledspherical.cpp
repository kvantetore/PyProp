#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>
#include <core/representation/combinedrepresentation.h>
#include <core/representation/coupledspherical/coupledsphericalharmonicrepresentation.h>

using namespace CoupledSpherical;

class CoupledSphericalSelectionRule
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

	CoupledSphericalSelectionRule() {}
	virtual ~CoupledSphericalSelectionRule() {}

	BasisPairList GetBasisPairs(const CoupledSphericalHarmonicRepresentation::Ptr angularRepresentation)
	{
		using namespace CoupledSpherical;

		//Setup angular index pair list
		typedef blitz::TinyVector<int, 2> BasisPair;
		typedef std::vector<BasisPair> BasisPairList;
		BasisPairList list;

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
		typedef BasisPairList::iterator Iterator;
		for (Iterator it=list.begin(); it!=list.end(); it++)
		{
		    indexPairs(pairIndex, 0) = (*it)(0);
		    indexPairs(pairIndex, 1) = (*it)(1);
		    pairIndex++;
		}


		return indexPairs;
	}

	virtual bool SelectionRule(CoupledSpherical::CoupledIndex const& left, CoupledSpherical::CoupledIndex const& right) = 0;
};


class CoupledSphericalSelectionRuleR12 : public CoupledSphericalSelectionRule
{
public:
	CoupledSpherical::ClebschGordan cg;

	CoupledSphericalSelectionRuleR12() {}
	virtual ~CoupledSphericalSelectionRuleR12() {}

	virtual bool SelectionRule(CoupledIndex const& left, CoupledIndex const& right)
	{
		int Lp = left.L;
		int Mp = left.M;
		int l1p = left.l1;
		int l2p = left.l2;
		int L = right.L;
		int M = right.M;
		int l1 = right.l1;
		int l2 = right.l2;

		bool nonzero = (L==Lp) && (M==Mp) && (L<=l1p+l2p) && (std::abs(l1p-l2p)<=L) && (L<=l1+l2) && (std::abs(l1-l2)<=L);

		return nonzero;
	}

};


class CoupledSphericalSelectionRuleLinearPolarizedField  : public CoupledSphericalSelectionRule
{
public:
	CoupledSphericalSelectionRuleLinearPolarizedField() {}
	virtual ~CoupledSphericalSelectionRuleLinearPolarizedField() {}

	virtual bool SelectionRule(CoupledIndex const& left, CoupledIndex const& right)
	{
		if ( (std::abs(left.L - right.L) == 1) && left.M == right.M)
		{
			if ( (std::abs(left.l1 - right.l1) == 1) && left.l2 == right.l2 )
			{
				return true;
			}

			if ( (std::abs(left.l2 - right.l2) == 1) && left.l1 == right.l1 )
			{
				return true;
			}
		}
		return false;
	}
};

class CoupledSphericalSelectionRuleDiagonal  : public CoupledSphericalSelectionRule
{
public:
	CoupledSphericalSelectionRuleDiagonal() {}
	virtual ~CoupledSphericalSelectionRuleDiagonal() {}

	virtual bool SelectionRule(CoupledIndex const& left, CoupledIndex const& right)
	{
		return (left.L == right.L) && (left.M == right.M) && (left.l1 == right.l1) && (left.l2 == right.l2);
	}
};



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


template<int Rank>
class CustomPotentialEvaluationR12 : public CustomPotentialCoupledSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

public:
	CustomPotentialEvaluationR12() {}
	virtual ~CustomPotentialEvaluationR12() {}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typedef CoupledSphericalHarmonicRepresentation CplHarmRepr;

		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CplHarmRepr::Ptr angRepr = boost::static_pointer_cast< CplHarmRepr >(repr->GetRepresentation(this->AngularRank));
	
		int r1Count = data.extent(this->RadialRank1);
		int r2Count = data.extent(this->RadialRank2);
		int angCount = data.extent(this->AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(this->RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(this->RadialRank2);

		ClebschGordan cg;

		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(this->RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		blitz::TinyVector<int, Rank> index;
		data = 0;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(this->AngularRank) = angIndex; 

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
				int lStop = std::max(std::max(l1, l1p), std::max(l2, l2p));
				for (int m1p=-lStop; m1p<=lStop; m1p++)
				{
					int m2p = M - m1p;
					for (int m1=-lStop; m1<=lStop; m1++)
					{
						int m2 = M - m1;
						int m3 = m1 - m1p;
						double cur = 1.0;
						cur *= cg(l1p, l2p, m1p, m2p, L, M);
						cur *= cg(l1, l2, m1, m2, L, M);
						cur *= cg(l1p, l1, -m1p, m1, l3, m3);
						cur *= cg(l2p, l2, -m2p, m2, l3, -m3);
						cur *= std::pow(-1., m1p + m2p + m3);
						l3Sum += cur;
					}
				}
				l3Sum *= l3Coeff; 

				for (int ri1=0; ri1<r1Count; ri1++)
				{
					double r1 = localr1(ri1);
					index(this->RadialRank1) = ri1;
			
					for (int ri2=0; ri2<r2Count; ri2++)
					{
						double r2 = localr2(ri2);
						index(this->RadialRank2) = ri2;
			
						double rmin = std::min(r1, r2);
						double rmax = std::max(r1, r2);
						double rfrac = rmin / rmax;
				
						data(index) += l3Sum * std::pow(rfrac, l3) / rmax;
					}
				}
			}
		}
	}
	
	static double CoefficientR12(int l1, int l2, int l3)
	{
		return 	std::sqrt( (2.*l1 + 1.) * (2.*l2 + 1.) / (2.*l3 + 1.) );
	}

	static double CondonShortleyPhase(int m)
	{
		if (m < 0) return 1.0;
		return std::pow(-1.0, m);
	}
};


/*
 * Potential evaluator for linearly polarized length gauge electric field
 */
template<int Rank>
class CustomPotentialEvaluationLinearPolarizedFieldLength : public CustomPotentialCoupledSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

public:
	CustomPotentialEvaluationLinearPolarizedFieldLength() {}
	virtual ~CustomPotentialEvaluationLinearPolarizedFieldLength() {}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typedef CoupledSphericalHarmonicRepresentation CplHarmRepr;

		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CplHarmRepr::Ptr angRepr = boost::static_pointer_cast< CplHarmRepr >(repr->GetRepresentation(this->AngularRank));
	
		int r1Count = data.extent(this->RadialRank1);
		int r2Count = data.extent(this->RadialRank2);
		int angCount = data.extent(this->AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(this->RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(this->RadialRank2);

		ClebschGordan cg;

		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(this->RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		blitz::TinyVector<int, Rank> index;
		data = 0;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(this->AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
			CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);
	
			if (std::abs(left.L - right.L) != 1) continue;
			if (left.M != right.M) continue;
	
			// "Left" quantum numbers
			int L = left.L;
			int M = left.M;
			int l1 = left.l1;
			int l2 = left.l2;
			
			// "Right" quantum numbers (Mp = M)
			int Lp = right.L;
			int l1p = right.l1;
			int l2p = right.l2;

			double I1 = 0;
			double I2 = 0;
			int lStop = std::max(std::max(l1, l1p), std::max(l2, l2p));

			for (int m1=-lStop; m1<=lStop; m1++)
			{
				int m2 = M - m1;
				int m1p = m1;
				int m2p = m2;

				if (std::abs(m1) > l1) continue;
				if (std::abs(m1p) > l1p) continue;
				if (std::abs(m2) > l2) continue;
				if (std::abs(m2p) > l2p) continue;

				double cur = cg(l1, l2, m1, m2, L, M) * cg(l1p, l2p, m1p, m2p, Lp, M);
				cur *= Coefficient(l1, l1p);
				cur *= cg(l1, l1p, 0, 0, 1, 0) * cg(l1, l1p, -m1, m1p, 1, 0);
				cur /= CondonShortleyPhase(m1);
				I1 += cur;
				
				cur = cg(l1, l2, m1, m2, L, M) * cg(l1p, l2p, m1p, m2p, Lp, M);
				cur *= Coefficient(l2, l2p);
				cur *= cg(l2, l2p, 0, 0, 1, 0) * cg(l2, l2p, -m2, m2p, 1, 0);
				cur /= CondonShortleyPhase(m2);
				I2 += cur;
			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				double r1 = localr1(ri1);
				index(this->RadialRank1) = ri1;

				for (int ri2=0; ri2<r2Count; ri2++)
				{
					double r2 = localr2(ri2);
					index(this->RadialRank2) = ri2;

					data(index) += (I1 * r1 + I2 * r2);
				}
			}
		}
	}

	static double Coefficient(int l, int lp)
	{
		return std::sqrt((2 * l + 1.0 ) * (2 * lp + 1.0) / (12 * M_PI));
	}

	static double CondonShortleyPhase(int m)
	{
		if (m < 0) return 1.0;
		return std::pow(-1.0, m);
	}
};


template<int Rank>
class CoupledSphericalKineticEnergyEvaluator : public CustomPotentialCoupledSphericalBase<Rank>
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

public:
	CoupledSphericalKineticEnergyEvaluator() {}
	virtual ~CoupledSphericalKineticEnergyEvaluator() {}

	double Mass;

	virtual void ApplyConfigSection(const ConfigSection &config)
	{
		CustomPotentialCoupledSphericalBase<Rank>::ApplyConfigSection(config);
		config.Get("mass", Mass);
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typedef CoupledSphericalHarmonicRepresentation CplHarmRepr;

		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CplHarmRepr::Ptr angRepr = boost::static_pointer_cast< CplHarmRepr >(repr->GetRepresentation(this->AngularRank));
	
		int r1Count = data.extent(this->RadialRank1);
		int r2Count = data.extent(this->RadialRank2);
		int angCount = data.extent(this->AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(this->RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(this->RadialRank2);

		ClebschGordan cg;

		BasisPairList angBasisPairs = GetBasisPairList(this->AngularRank);

		if (data.extent(this->RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(this->RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(this->AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		blitz::TinyVector<int, Rank> index;
		data = 0;
	
		for (int ri1=0; ri1<r1Count; ri1++)
		{
			index(this->RadialRank1) = ri1;
			double r1 = localr1(ri1);
	
			for (int ri2=0; ri2<r2Count; ri2++)
			{
				index(this->RadialRank2) = ri2;
				double r2 = localr2(ri2);
	
				for (int angIndex=0; angIndex<angCount; angIndex++)
				{
					index(this->AngularRank) = angIndex;

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
	
					data(index) = V1 + V2;
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
	int RadialRank1;
	int RadialRank2;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("z", Z);
		config.Get("radial_rank1", RadialRank1);
		config.Get("radial_rank2", RadialRank2);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r1 = pos(RadialRank1);
		double r2 = pos(RadialRank2);
		return - Z/r1 - Z/r2;
	}
};



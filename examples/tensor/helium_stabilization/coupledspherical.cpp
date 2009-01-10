#include "coupledbase.h"

/*
 * Repulsion Energy 1/|r1 - r2|
 */
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
 * Angular Kinetic Energy
 */
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




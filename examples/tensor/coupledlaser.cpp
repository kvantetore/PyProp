#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>
#include <core/representation/combinedrepresentation.h>
#include <core/representation/coupledspherical/coupledsphericalharmonicrepresentation.h>

#include "laserhelper.h"

using namespace CoupledSpherical;

/* First part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics
 *
 * <Ylm | - \frac{1}{r} \sin \theta \partialdiff{}{\theta} 
 *	      - \frac{\cos \theta}{r} | Yl'm'>
 */
template<int Rank>
class CustomPotential_LaserVelocity_CoupledSpherical
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank1;
	int RadialRank2;

public:
	CustomPotential_LaserVelocity_CoupledSpherical() {}
	virtual ~CustomPotential_LaserVelocity_CoupledSpherical() {}

	void ApplyConfigSection(const ConfigSection &config)
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

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CoupledSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< CoupledSphericalHarmonicRepresentation >(repr->GetRepresentation(AngularRank));
	
		int r1Count = data.extent(RadialRank1);
		int r2Count = data.extent(RadialRank2);
		int angCount = data.extent(AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(RadialRank2);
		BasisPairList angBasisPairs = GetBasisPairList(AngularRank);

		ClebschGordan cg;

		if (data.extent(RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		cplx IM(0,1.0);

		data = 0;
		blitz::TinyVector<int, Rank> index;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
			CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);

			//"Left" quantum numbers
			int l1 = left.l1;
			int l2 = left.l2;
			int L = left.L;
			int M = left.M;
			
			//"Right" quantum numbers (Mp = M)
			int l1p = right.l1;
			int l2p = right.l2;
			int Lp = right.L;
			int Mp = right.M;

			//Rough selection rule
			if (M != Mp) continue;
			int lStop = std::max(std::max(l1, l1p), std::max(l2, l2p));

			double coupling1 = 0;
			double coupling2 = 0;

			for (int m=-lStop; m<=lStop; m++)
			{
				double curCg = cg(l1, l2, m, M-m, L, M) * cg(l1p, l2p, m, Mp-m, Lp, Mp);

				if (LaserHelper::kronecker(l2, l2p) != 0)
				{
					if (std::abs(m) <= l1 && std::abs(m) <= l1p)
					{
						double C = LaserHelper::C(l1p, m) * LaserHelper::kronecker(l1, l1p+1);
						double D = LaserHelper::D(l1p, m) * LaserHelper::kronecker(l1, l1p-1);
						double E = LaserHelper::E(l1p, m) * LaserHelper::kronecker(l1, l1p+1);
						double F = LaserHelper::F(l1p, m) * LaserHelper::kronecker(l1, l1p-1);

						coupling1 += (- (C + D) - (E + F)) * curCg;
					}
				}

				if (LaserHelper::kronecker(l1, l1p) != 0)
				{
					if (std::abs(m) <= l2 && std::abs(m) <= l2p)
					{
						double C = LaserHelper::C(l2p, M-m) * LaserHelper::kronecker(l2, l2p+1);
						double D = LaserHelper::D(l2p, M-m) * LaserHelper::kronecker(l2, l2p-1);
						double E = LaserHelper::E(l2p, M-m) * LaserHelper::kronecker(l2, l2p+1);
						double F = LaserHelper::F(l2p, M-m) * LaserHelper::kronecker(l2, l2p-1);

						coupling2 += (- (C + D) - (E + F)) * curCg;
					}
				}
			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				index(RadialRank1) = ri1;
				double r1 = localr1(ri1);

				for (int ri2=0; ri2<r2Count; ri2++)
				{
					index(RadialRank2) = ri2;
					double r2 = localr2(ri2);
				
					data(index) = - IM * (coupling1/r1 + coupling2/r2);
				}
			}
		}
	}
};

/* Second part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics.
 *
 * Should be used with first order differentiation in r1
 *
 */
template<int Rank>
class CustomPotential_LaserVelocityDerivativeR1_CoupledSpherical
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank1;
	int RadialRank2;

public:
	CustomPotential_LaserVelocityDerivativeR1_CoupledSpherical() {}
	virtual ~CustomPotential_LaserVelocityDerivativeR1_CoupledSpherical() {}

	void ApplyConfigSection(const ConfigSection &config)
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

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CoupledSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< CoupledSphericalHarmonicRepresentation >(repr->GetRepresentation(AngularRank));
	
		int r1Count = data.extent(RadialRank1);
		int r2Count = data.extent(RadialRank2);
		int angCount = data.extent(AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(RadialRank2);
		BasisPairList angBasisPairs = GetBasisPairList(AngularRank);

		ClebschGordan cg;

		if (data.extent(RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		cplx IM(0,1.0);

		data = 0;
		blitz::TinyVector<int, Rank> index;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
			CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);

			//"Left" quantum numbers
			int l1 = left.l1;
			int l2 = left.l2;
			int L = left.L;
			int M = left.M;
			
			//"Right" quantum numbers (Mp = M)
			int l1p = right.l1;
			int l2p = right.l2;
			int Lp = right.L;
			int Mp = right.M;

			//Rough selection rule
			if (M != Mp) continue;
			int lStop = std::max(std::max(l1, l1p), std::max(l2, l2p));

			double coupling1 = 0;
			double coupling2 = 0;

			for (int m=-lStop; m<=lStop; m++)
			{
				double curCg = cg(l1, l2, m, M-m, L, M) * cg(l1p, l2p, m, Mp-m, Lp, Mp);

				if (LaserHelper::kronecker(l2, l2p) != 0)
				{
					if (std::abs(m) <= l1 && std::abs(m) <= l1p)
					{
						double E = LaserHelper::E(l1p, m) * LaserHelper::kronecker(l1, l1p+1);
						double F = LaserHelper::F(l1p, m) * LaserHelper::kronecker(l1, l1p-1);

						coupling1 += (E + F) * curCg;
					}
				}
			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				index(RadialRank1) = ri1;
				double r1 = localr1(ri1);

				for (int ri2=0; ri2<r2Count; ri2++)
				{
					index(RadialRank2) = ri2;
					double r2 = localr2(ri2);
				
					data(index) = - IM * (coupling1);
				}
			}
		}
	}
};

/* Third part of the linearly polarized laser in the velocity gauge
 * expressed in spherical harmonics.
 *
 * Should be used with first order differentiation in r2
 *
 */
template<int Rank>
class CustomPotential_LaserVelocityDerivativeR2_CoupledSpherical
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;
	int AngularRank;
	int RadialRank1;
	int RadialRank2;

public:
	CustomPotential_LaserVelocityDerivativeR2_CoupledSpherical() {}
	virtual ~CustomPotential_LaserVelocityDerivativeR2_CoupledSpherical() {}

	void ApplyConfigSection(const ConfigSection &config)
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

	BasisPairList GetBasisPairList(int rank)
	{
		if (rank == AngularRank)
			return AngularBasisPairs;
		else
			return BasisPairList();
	}

	virtual void UpdatePotentialData(typename blitz::Array<cplx, Rank> data, typename Wavefunction<Rank>::Ptr psi, cplx timeStep, double curTime)
	{
		using namespace CoupledSpherical;

		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast< CmbRepr >(psi->GetRepresentation());
		CoupledSphericalHarmonicRepresentation::Ptr angRepr = boost::static_pointer_cast< CoupledSphericalHarmonicRepresentation >(repr->GetRepresentation(AngularRank));
	
		int r1Count = data.extent(RadialRank1);
		int r2Count = data.extent(RadialRank2);
		int angCount = data.extent(AngularRank);

		blitz::Array<double, 1> localr1 = psi->GetRepresentation()->GetLocalGrid(RadialRank1);
		blitz::Array<double, 1> localr2 = psi->GetRepresentation()->GetLocalGrid(RadialRank2);
		BasisPairList angBasisPairs = GetBasisPairList(AngularRank);

		ClebschGordan cg;

		if (data.extent(RadialRank1) != r1Count) throw std::runtime_error("Invalid r1 size");
		if (data.extent(RadialRank2) != r2Count) throw std::runtime_error("Invalid r2 size");
		if (data.extent(AngularRank) != angBasisPairs.extent(0)) throw std::runtime_error("Invalid ang size");

		cplx IM(0,1.0);

		data = 0;
		blitz::TinyVector<int, Rank> index;
	
		for (int angIndex=0; angIndex<angCount; angIndex++)
		{
			index(AngularRank) = angIndex;

			int leftIndex = angBasisPairs(angIndex, 0);
			int rightIndex = angBasisPairs(angIndex, 1);
	
			CoupledIndex left = angRepr->Range.GetCoupledIndex(leftIndex);
			CoupledIndex right = angRepr->Range.GetCoupledIndex(rightIndex);

			//"Left" quantum numbers
			int l1 = left.l1;
			int l2 = left.l2;
			int L = left.L;
			int M = left.M;
			
			//"Right" quantum numbers (Mp = M)
			int l1p = right.l1;
			int l2p = right.l2;
			int Lp = right.L;
			int Mp = right.M;

			//Rough selection rule
			if (M != Mp) continue;
			int lStop = std::max(std::max(l1, l1p), std::max(l2, l2p));

			double coupling2 = 0;

			for (int m=-lStop; m<=lStop; m++)
			{
				double curCg = cg(l1, l2, m, M-m, L, M) * cg(l1p, l2p, m, Mp-m, Lp, Mp);

				if (LaserHelper::kronecker(l1, l1p) != 0)
				{
					if (std::abs(m) <= l2 && std::abs(m) <= l2p)
					{
						double E = LaserHelper::E(l2p, m) * LaserHelper::kronecker(l2, l2p+1);
						double F = LaserHelper::F(l2p, m) * LaserHelper::kronecker(l2, l2p-1);

						coupling2 += (E + F) * curCg;
					}
				}
			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				index(RadialRank1) = ri1;
				double r1 = localr1(ri1);

				for (int ri2=0; ri2<r2Count; ri2++)
				{
					index(RadialRank2) = ri2;
					double r2 = localr2(ri2);
				
					data(index) = - IM * (coupling2);
				}
			}
		}
	}
};

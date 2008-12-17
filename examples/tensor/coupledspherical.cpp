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
		bool nonzero = left.L == right.L && left.M == right.M;
		if (nonzero)
		{
			nonzero = false;

			int L = left.L;
			int M = left.M;
			int l1p = left.l1;
			int l2p = left.l2;
			int l1 = right.l1;
			int l2 = right.l2;
		
			int minL3 = std::max(std::abs(l1 - l1p), std::abs(l2 - l2p));
			int maxL3 = std::min(l1+l1p, l2+l2p);

			double eps = 1e-12;
		
			for (int l3=minL3; l3<=maxL3; l3++)
			{
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
				
				nonzero = nonzero || std::abs(l3Sum)>eps;
			}
		}

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
		if ( (std::abs(left.L - right.L) == 1) && left.M == right.M )
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


template<int Rank>
class CustomPotentialEvaluationR12
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;

public:
	CustomPotentialEvaluationR12() {}
	virtual ~CustomPotentialEvaluationR12() {}

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

/*
 * Helper functions for coupled spherical harmonics laser interaction matrix elements
 */
class HelperFunctions
{
	public:

	static double CoefficientC(int l, int m)
	{
		//Avoid division by zero or negative square root
		if (l - 1 <= 0) return 0.0;

		double C = std::sqrt( (2.0 * l + 1.0) / (2.0 * l - 1.0) * (l - m) * (l + m) );
		//cout << "l = " << l << " m = " << m << " C = " << C << endl;
		return C;
	}

	static double CoefficientD(int l, int m)
	{
		double D =  (l - m + 1.0) / CondonShortleyPhase(m);
		D *= std::sqrt( (4 * M_PI) / (2.0 * l + 3.0) * StableFactorial(l + m + 1, l - m + 1) );
		return D;
	}

	static double CoefficientE(int l, int m)
	{
		//Avoid division by zero or negative square root
		if (2 * l - 1 <= 0) return 0.0;

		double E = (l + m) / CondonShortleyPhase(m);
		E *= std::sqrt( (4 * M_PI) / (2.0 * l - 1.0) * StableFactorial(l + m - 1, l - m - 1) );
		return E;
	}

	static double CondonShortleyPhase(int m)
	{
		return std::pow(-1.0, m);
	}

	static double Kronecker(int a, int b)
	{
		if (a == b) return 1.0;
		return 0.0;
	}

	static double StableFactorial(int numerator, int denominator)
	{
		double numeratorLogFactorial = 0;
		for (int a = numerator; a>0; a--)
		{
			numeratorLogFactorial += std::log(a);
		}

		double denominatorLogFactorial = 0;
		for (int b = numerator; b>0; b--)
		{
			denominatorLogFactorial += std::log(b);
		}

		return std::exp(numeratorLogFactorial - denominatorLogFactorial);
	}
};

/*
 * Potential evaluator for linearly polarized velocity gauge electric field,
 * 1/r1 and 1/r2 term.
 */
template<int Rank>
class CustomPotentialEvaluationLinearPolarizedField
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;

public:
	CustomPotentialEvaluationLinearPolarizedField() {}
	virtual ~CustomPotentialEvaluationLinearPolarizedField() {}

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
	
			if (std::abs(left.L - right.L) != 1) continue;
			if (left.M != right.M) continue;
	
			int L = left.L;
			int M = left.M;
			int l1 = left.l1;
			int l2 = left.l2;

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

				double cur = HelperFunctions::CoefficientD(l1, m1) * HelperFunctions::Kronecker(l1, l1p + 1);
				cur += HelperFunctions::CoefficientE(l1, m1) * HelperFunctions::Kronecker(l1, l1p - 1);
				cur *= l1p;
				cur -= HelperFunctions::CoefficientC(l1, m1) * HelperFunctions::Kronecker(l1, l1p - 1);
				cur *= cg(l1p, l2p, m1p, m2p, Lp, M);
				cur *= cg(l1, l2, m1, m2, L, M);
				cur *= HelperFunctions::Kronecker(l2, l2p);
				I1 += cur;

				cur = HelperFunctions::CoefficientD(l2, m2) * HelperFunctions::Kronecker(l2, l2p + 1);
				cur += HelperFunctions::CoefficientE(l2, m2) * HelperFunctions::Kronecker(l2, l2p - 1);
				cur *= l2p;
				cur -= HelperFunctions::CoefficientC(l2, m2) * HelperFunctions::Kronecker(l2, l2p - 1);
				cur *= cg(l1p, l2p, m1p, m2p, Lp, M);
				cur *= cg(l1, l2, m1, m2, L, M);
				cur *= HelperFunctions::Kronecker(l1, l1p);
				I2 += cur;
			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				double r1 = localr1(ri1);

				for (int ri2=0; ri2<r2Count; ri2++)
				{
					double r2 = localr2(ri2);

					data(ri1, ri2, angIndex) += cplx(0., .1) * (I1 / r1 + I2 / r2);
				}
			}
		}
	}
};

/*
 * Potential evaluator for linearly polarized velocity gauge electric field,
 * radial derivative in direction 1.
 */
template<int Rank>
class CustomPotentialEvaluationLinearPolarizedFieldDerivativeR1
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;

public:
	CustomPotentialEvaluationLinearPolarizedFieldDerivativeR1() {}
	virtual ~CustomPotentialEvaluationLinearPolarizedFieldDerivativeR1() {}

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
	
			if (std::abs(left.L - right.L) != 1) continue;
			if (left.M != right.M) continue;
	
			int L = left.L;
			int M = left.M;
			int Lp = right.L;
			int l1p = left.l1;
			int l2p = left.l2;
			int l1 = right.l1;
			int l2 = right.l2;

			int lStop = std::max(std::max(l1, l1p), std::max(l2, l2p));

			double I1 = 0;
			for (int m1=-lStop; m1<=lStop; m1++)
			{
				int m2 = M - m1;
				int m1p = m1;
				int m2p = m2;

				if (std::abs(m1) > l1) continue;
				if (std::abs(m1p) > l1p) continue;
				if (std::abs(m2) > l2) continue;
				if (std::abs(m2p) > l2p) continue;

				double cur = HelperFunctions::CoefficientD(l1, m1) * HelperFunctions::Kronecker(l1, l1p+1);
				cur += HelperFunctions::CoefficientE(l1, m1) * HelperFunctions::Kronecker(l1, l1p-1);
				cur *= cg(l1p, l2p, m1p, m2p, Lp, M);
				cur *= cg(l1, l2, m1, m2, L, M);
				cur *= HelperFunctions::Kronecker(l2, l2p);
				I1 += cur;
			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				for (int ri2=0; ri2<r2Count; ri2++)
				{
					data(ri1, ri2, angIndex) -= cplx(0., .1) * I1;
				}
			}
		}
	}
};

/*
 * Potential evaluator for linearly polarized velocity gauge electric field,
 * radial derivative in direction 2.
 */

template<int Rank>
class CustomPotentialEvaluationLinearPolarizedFieldDerivativeR2
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

private:
	BasisPairList AngularBasisPairs;

public:
	CustomPotentialEvaluationLinearPolarizedFieldDerivativeR2() {}
	virtual ~CustomPotentialEvaluationLinearPolarizedFieldDerivativeR2() {}

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
	
			if (std::abs(left.L - right.L) != 1) continue;
			if (left.M != right.M) continue;
	
			int L = left.L;
			int M = left.M;
			int Lp = right.L;
			int l1p = left.l1;
			int l2p = left.l2;
			int l1 = right.l1;
			int l2 = right.l2;

			int lStop = std::max(std::max(l1, l1p), std::max(l2, l2p));

			double I2 = 0;
			for (int m1=-lStop; m1<=lStop; m1++)
			{
				int m2 = M - m1;
				int m1p = m1;
				int m2p = m2;

				if (std::abs(m1) > l1) continue;
				if (std::abs(m1p) > l1p) continue;
				if (std::abs(m2) > l2) continue;
				if (std::abs(m2p) > l2p) continue;

				double cur = HelperFunctions::CoefficientD(l2, m2) * HelperFunctions::Kronecker(l2, l2p+1);
				cur += HelperFunctions::CoefficientE(l2, m2) * HelperFunctions::Kronecker(l2, l2p-1);
				cur *= cg(l1p, l2p, m1p, m2p, Lp, M);
				cur *= cg(l1, l2, m1, m2, L, M);
				cur *= HelperFunctions::Kronecker(l1, l1p);
				I2 += cur;
			}

			for (int ri1=0; ri1<r1Count; ri1++)
			{
				for (int ri2=0; ri2<r2Count; ri2++)
				{
					data(ri1, ri2, angIndex) -= cplx(0., .1) * I2;
				}
			}
		}
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
	virtual ~CoupledSphericalKineticEnergyEvaluator() {}

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



#ifndef COUPLEDSPHERICALSELECTIONRULE_H
#define COUPLEDSPHERICALSELECTIONRULE_H

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


class CoupledSphericalSelectionRuleR12Old : public CoupledSphericalSelectionRule
{

public:
	CoupledSpherical::ClebschGordan cg;

	CoupledSphericalSelectionRuleR12Old() {}

	virtual ~CoupledSphericalSelectionRuleR12Old() {}

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
		

class CoupledSphericalSelectionRuleR12 : public CoupledSphericalSelectionRule
{
private:
	int MultipoleCutoff;

public:
	CoupledSpherical::ClebschGordan cg;

	CoupledSphericalSelectionRuleR12(): MultipoleCutoff(std::numeric_limits<int>::max()) {}

	CoupledSphericalSelectionRuleR12(int multipoleCutoff)
	{
		MultipoleCutoff = multipoleCutoff;
		std::cout << "Multipole cutoff = " << MultipoleCutoff << std::endl;
	}

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
		
		int minL3 = std::max(std::abs(l1 - l1p), std::abs(l2 - l2p));
		int maxL3 = std::min(l1+l1p, l2+l2p);
		maxL3 = std::min(MultipoleCutoff, maxL3);

		double l3Sum = 0.0;
		for (int l3=minL3; l3<=maxL3; l3++)
		{
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
		}

		return nonzero && (std::abs(l3Sum) > 1e-14);
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


#endif


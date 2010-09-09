#ifndef SPHERICALBASISSELECTIONRULE_H
#define SPHERICALBASISSELECTIONRULE_H

#include <core/representation/sphericalbasis/sphericalharmonicbasisrepresentation.h>

using namespace SphericalBasis;

class SphericalBasisSelectionRule
{
public:
	typedef blitz::Array<int, 2> BasisPairList;

	SphericalBasisSelectionRule() {}
	virtual ~SphericalBasisSelectionRule() {}

	BasisPairList GetBasisPairs(const SphericalHarmonicBasisRepresentation::Ptr angularRepresentation)
	{
		using namespace SphericalBasis;

		//Setup angular index pair list
		typedef blitz::TinyVector<int, 2> BasisPair;
		typedef std::vector<BasisPair> BasisPairList;
		BasisPairList list;

		int angCount = angularRepresentation->Range.Count();
		for (int leftIndex=0; leftIndex<angCount; leftIndex++)
		{
			LmIndex left = angularRepresentation->Range.GetLmIndex(leftIndex);
			
			for (int rightIndex=0; rightIndex<angCount; rightIndex++)
			{
			LmIndex right = angularRepresentation->Range.GetLmIndex(rightIndex);
				
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

	virtual bool SelectionRule(SphericalBasis::LmIndex const& left, SphericalBasis::LmIndex const& right) = 0;
};


/*
 * Selection rule for linearly polarized field in the z direction
 */
class SphericalBasisSelectionRuleLinearPolarizedField  : public SphericalBasisSelectionRule
{
public:
	SphericalBasisSelectionRuleLinearPolarizedField() {}
	virtual ~SphericalBasisSelectionRuleLinearPolarizedField() {}

	virtual bool SelectionRule(LmIndex const& left, LmIndex const& right)
	{
		if ( (std::abs(left.l - right.l) == 1) && left.m == right.m )
		{
			return true;
		}
		return false;
	}
};

/*
 * Diagonal selection rule (l == l', m == m')
 */
class SphericalBasisSelectionRuleDiagonal  : public SphericalBasisSelectionRule
{
public:
	SphericalBasisSelectionRuleDiagonal() {}
	virtual ~SphericalBasisSelectionRuleDiagonal() {}

	virtual bool SelectionRule(LmIndex const& left, LmIndex const& right)
	{
		return (left.l == right.l) && (left.m == right.m);
	}
};

/*
 * Selection rule for linearly polarized field perpendicular to the z-axis (x- and y-fields)
 */
class SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular  : public SphericalBasisSelectionRule
{
public:
	SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular() {}
	virtual ~SphericalBasisSelectionRuleLinearPolarizedFieldPerpendicular() {}

	virtual bool SelectionRule(LmIndex const& left, LmIndex const& right)
	{
		if ( (std::abs(left.l - right.l) == 1) && (std::abs(left.m - right.m) == 1))
		{
			return true;
		}
		return false;
	}
};

/*
 * Selection rule for linearly polarized field perpendicular at an angle to the z-axis.
 */
class SphericalBasisSelectionRuleLinearPolarizedFieldAngle  : public SphericalBasisSelectionRule
{
public:
	SphericalBasisSelectionRuleLinearPolarizedFieldAngle() {}
	virtual ~SphericalBasisSelectionRuleLinearPolarizedFieldAngle() {}

	virtual bool SelectionRule(LmIndex const& left, LmIndex const& right)
	{
		if ( (std::abs(left.l - right.l) == 1) && (std::abs(left.m - right.m) <= 1))
		{
			return true;
		}
		return false;
	}
};

#endif


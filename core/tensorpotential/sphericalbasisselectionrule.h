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


//
// Selection rule for diatomic Coulomb potentials
//
class SphericalBasisSelectionRuleDiatomicCoulomb : public SphericalBasisSelectionRule
{
	private:
		//a max limit for max l that should correspond to number of terms in
		//multipole expansion
		int MultipoleCutoff;  
	
	public:
		//ClebshGordan Coefficent calculator
		ClebschGordan cg;
	
		//constructors
		SphericalBasisSelectionRuleDiatomicCoulomb() :
			 MultipoleCutoff(std::numeric_limits<int>::max()) {}
	
		SphericalBasisSelectionRuleDiatomicCoulomb(int multipoleCutoff)
		{
			MultipoleCutoff = multipoleCutoff;
			#ifdef PYPROP_DEBUG
			std::cout << "Multipole cutoff = " << MultipoleCutoff << std::endl;
			#endif
		}

		//destructor
		virtual ~SphericalBasisSelectionRuleDiatomicCoulomb() {}
		
		//Method that is called for two coupled indicies left anf rigth
		//returns true if coupling is non-zero
		virtual bool SelectionRule(LmIndex const & left,
			 LmIndex const & right)
		{
			int mp = left.m;
			int lp = left.l;
			
			int m = right.m;
			int l = right.l;

			//checks that the given coupling indices are legal
			int minL3 = std::abs(l - lp);
			int maxL3 = l + lp;
			maxL3 = std::min(MultipoleCutoff, maxL3);

			//Sum of all contributions to the coupling
			double l3Sumfinal = 0;

			for(int l3 = minL3; l3<=maxL3; l3++)
			{
				if(l3 % 2 == 1) continue;

				double l3Coeff = 1.0;
				l3Coeff *= Coefficient(l,lp);
				l3Coeff *= cg(l, l3, 0, 0, lp, 0);
			
				double l3Sum = 0;
				
				for(int m3 = -l3; m3 <= l3; m3++)
				{
					double cur = 1; 

					cur *= 2.0;
					cur *=CondonShortleyPhase(-m3);
					cur *= MultipoleCoeff(l3);
						
					cur *= cg(l,l3,m,m3,lp,mp);
						
					l3Sum += cur;
				}
				l3Sumfinal += l3Sum*l3Coeff;
			}
		return (std::abs(l3Sumfinal) > 1e-14);
	}

	static double Coefficient(int a, int b)
	{
		return std::sqrt((2. *a + 1.) / (2. *b +1.));
	}

	static double MultipoleCoeff(int c)
	{
		return std::sqrt((4. * M_PI) / (2. * c + 1.));
	}
	
	static double CondonShortleyPhase(int m)
	{
		if(m < 0) return 1.0;
		return std::pow(-1.0, m);
	}
	
	static int kronecker(int a, int b)
	{
		if(a == b)
			return 1;
		else
			return 0;
	}
};

#endif


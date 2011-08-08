#ifndef COUPLEDSPHERICALSELECTIONRULE_H
#define COUPLEDSPHERICALSELECTIONRULE_H

#include "src/representation/coupledsphericalharmonicrepresentation.h"

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
		#ifdef PYPROP_DEBUG
		std::cout << "Multipole cutoff = " << MultipoleCutoff << std::endl;
		#endif
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
			double l3SumInner = 0.0;
			double l3Coeff = cg(l1p, l1, 0, 0, l3, 0);
			l3Coeff *= cg(l2p, l2, 0, 0, l3, 0);
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
					l3SumInner += cur;
				}
			}
			l3Sum += l3SumInner * l3Coeff;
			//l3Sum += l3SumInner;
		}

		return nonzero && (std::abs(l3Sum) > 1e-14);
	}
};

//Added ----------------------------------------------------------------
class CoupledSphericalSelectionRuleDiatomicCoulomb : public 
		CoupledSphericalSelectionRule
{
	private:
		//a max limit for max L that should correspond to number of terms in
		//multipole expansion
		int MultipoleCutoff;  
	
	public:
		//ClebshGordan Coeficant calculator
		CoupledSpherical::ClebschGordan cg;
	
		//constructors
		CoupledSphericalSelectionRuleDiatomicCoulomb() :
			 MultipoleCutoff(std::numeric_limits<int>::max()) {}
	
		CoupledSphericalSelectionRuleDiatomicCoulomb(int multipoleCutoff)
		{
			MultipoleCutoff = multipoleCutoff;
			#ifdef PYPROP_DEBUG
			std::cout << "Multipole cutoff = " << MultipoleCutoff << std::endl;
			#endif
		}

		//destructor
		virtual ~CoupledSphericalSelectionRuleDiatomicCoulomb() {}
		
		//Method that is called for two coupled indicies left anf rigth
		//returns true if coupling is non-zero
		virtual bool SelectionRule(CoupledIndex const & left,
			 CoupledIndex const & right)
		{
			int Lp = left.L;
			int Mp = left.M;
			int l1p = left.l1;
			int l2p = left.l2;
			
			int L = right.L;
			int M = right.M;
			int l1 = right.l1;
			int l2 = right.l2;

			//checks that the given coupling indices are legal
			//bool nonzero = (L==Lp) && (M == Mp) && (Lp <= l1p + l2p) &&
			//	(std::abs(l1p-l2p) <= Lp) && (L <= l1 +l2) && (std::abs(l1-l2)
			//	<= L);
				
			int minL3 = std::abs(l1 - l1p);
			int maxL3 = l1 + l1p;
			maxL3 = std::min(MultipoleCutoff, maxL3);

			//Sum of all contributions to the coupling
			double l3Sumfinal = 0;

			for(int l3 = minL3; l3<=maxL3; l3++)
			{
				double l3Coeff = 1.0;
				l3Coeff *= Coefficient(l1,l1p);
				l3Coeff *= cg(l1, l3, 0, 0, l1p, 0);
				l3Coeff *= kronecker(l2,l2p);
			
				double l3Sum = 0;
				
				for(int m1p = -l1p; m1p <= l1p; m1p++)
				{
					int m2p = Mp - m1p;
					for(int m1 = -l1; m1 <= l1; m1++)
					{
						int m2 = M - m1;
						int m3 = m1p - m1;
						
						if(l3 % 2 ==1) continue;
						if(std::abs(m1) > l1) continue;
						if(std::abs(m1p)>l1p) continue;
						if(std::abs(m2) > l2) continue;
						if(std::abs(m2p)>l2p) continue;
						if(std::abs(m3) > l3) continue;
			
						double cur = 1; 

						cur *= 2.0;
						cur *=CondonShortleyPhase(-m3);
						cur *= MultipoleCoeff(l3);
						
						cur *= cg(l1p,l2p,
							m1p,m2p,Lp,Mp);
						cur *= cg(l1,l2,m1,m2,L,M);
						cur *= cg(l1,l3,m1,m3,l1p,m1p);
						
						cur *= kronecker
							(m2, m2p);
						
						l3Sum += cur;
					}
				}
				l3Sumfinal += l3Sum*l3Coeff;
			}
				
			// r_1 <=> r_2
	
			minL3 = std::abs(l2 - l2p);
			maxL3 = l2 + l2p;
			maxL3 = std::min(MultipoleCutoff,maxL3);		

			for(int l3 = minL3; l3 <= maxL3; l3++)
			{
				double l3Coeff = 1.0;
				l3Coeff *= Coefficient(l2,l2p);
				l3Coeff *= cg(l2,l3,0,0,l2p,0);
				l3Coeff *= kronecker(l1,l1p);
			
				double l3Sum = 0;
				
				for(int m1p = -l1p; m1p <= l1p; m1p++)
				{
					int m2p = Mp - m1p;
					
					for(int m1= -l1; m1 <= l1; m1++)
					{
						int m2 = M - m1;
						int m3 = m2p - m2;

						//makes sure that only even ls are non-zero
						if((l3) % 2 == 1) continue;
						if(std::abs(m1)>l1) continue;
						if(std::abs(m1p)>l1p) continue;
						if(std::abs(m2) > l2) continue;
						if(std::abs(m2p)>l2p) continue;
						if(std::abs(m3) > l3) continue;
					
						
						double cur = 1; 
						cur *= CondonShortleyPhase(-m3);
						cur *= MultipoleCoeff(l3);
						cur *= 2.0;
						cur *=cg(l1p,l2p,m1p,m2p,Lp,Mp);
						cur *= cg(l1,l2,m1,m2,L,M);
						cur *= cg(l2,l3,m2,m3,l2p,m2p);
						
						cur *= kronecker
							(m1,m1p);
					
						l3Sum += cur;
					}
				}
			
				l3Sumfinal += l3Sum *= l3Coeff;
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

//adding stop ----------------------------------------------------------------
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

/*
 * Selection rule for linearly polarized field perpendicular to the z-axis (x- and y-fields)
 */
class CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular  : public CoupledSphericalSelectionRule
{
public:
	CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular() {}
	virtual ~CoupledSphericalSelectionRuleLinearPolarizedFieldPerpendicular() {}

	virtual bool SelectionRule(CoupledIndex const& left, CoupledIndex const& right)
	{
		if ( (std::abs(left.L - right.L) == 1) && (std::abs(left.M - right.M) == 1))
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

/*
 * Selection rule for linearly polarized field perpendicular at an angle to the z-axis.
 */
class CoupledSphericalSelectionRuleLinearPolarizedFieldAngle  : public CoupledSphericalSelectionRule
{
public:
	CoupledSphericalSelectionRuleLinearPolarizedFieldAngle() {}
	virtual ~CoupledSphericalSelectionRuleLinearPolarizedFieldAngle() {}

	virtual bool SelectionRule(CoupledIndex const& left, CoupledIndex const& right)
	{
		if ( (std::abs(left.L - right.L) == 1) && (std::abs(left.M - right.M) <= 1))
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

#endif


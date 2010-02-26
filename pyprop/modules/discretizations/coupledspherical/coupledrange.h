#ifndef COUPLEDRANGE_H
#define COUPLEDRANGE_H

//#include <boost/unordered_map.hpp>
#include <map>
//#include <boost/hash.hpp>
//#include <boost/foreach.hpp>

#include "../../common.h"
#include "../orthogonalrepresentation.h"
#include <gsl/gsl_sf_coupling.h>


namespace CoupledSpherical
{

class ClebschGordan
{
public:
	ClebschGordan() {}

	double operator()(int l1, int l2, int m1, int m2, int L, int M)
	{
		gsl_sf_result result;
		int status = gsl_sf_coupling_3j_e(2*l1, 2*l2, 2*L, 2*m1, 2*m2, -2*M, &result);
		if (status != 0)
		{
			cout << "Error finding cg coeff " << l1 << ", " << l2 << ", " << m1 << ", " << m2 << ", " << L << ", " << M << endl;
		}
		return pow(-1., M + l1 - l2) * std::sqrt(2*L + 1.0) * result.val;
	}
};

class CoupledIndex
{
public:
	typedef int IndexType;

	IndexType l1;
	IndexType l2;
	IndexType L;
	IndexType M;

	CoupledIndex() : l1(-1), l2(-1), L(-1), M(-1) {}
	CoupledIndex(IndexType l1, IndexType l2, IndexType L, IndexType M) : l1(l1), l2(l2), L(L), M(M) {}

	bool operator==(CoupledIndex const& other) const
	{
		return l1 == other.l1 
			&& l2 == other.l2
			&& L == other.L
			&& M == other.M;
	}

	unsigned int GetInt() const
	{
		typedef unsigned int u;
		return ((u)l1<<7*3) + ((u)l2<<7*2) + ((u)L<<7) + M;
	}
	
};

/*
inline std::size_t hash_value(CoupledIndex const& p)
{
	std::size_t seed = 0;
	boost::hash_combine(seed, p.l1);
	boost::hash_combine(seed, p.l2);
	boost::hash_combine(seed, p.L);
	boost::hash_combine(seed, p.M);
	return seed;
}
*/

class CoupledIndexCompareLT
{
public:
	bool operator()(CoupledIndex const& p1, CoupledIndex const& p2) const
	{
		return p1.GetInt() < p2.GetInt();
	}
};

//typedef boost::unordered_map<CoupledIndex, int> CoupledIndexMap;
typedef std::map<CoupledIndex, int, CoupledIndexCompareLT> CoupledIndexMap;

class CoupledRange
{
private:
	blitz::Array<double, 1> Grid;		//Grid is a double representation of the index, it is only 
										//used to satisfy the Representation interface
	blitz::Array<double, 1> Weights;    //Weights is the Clebsch-Gordan coefficients 

	blitz::Array<CoupledIndex, 1> IndexList;  //Maps from a grid index to a coupled index
	CoupledIndexMap IndexMap;				  //Maps from a coupled index to a grid index

	int curGridIndex;
	
public:
	//Constructors
	CoupledRange() : curGridIndex(0) {}
	~CoupledRange() {}

	//Returns the number possible lm values
	inline int Count() const
	{
		return Grid.extent(0);	
	}

	const blitz::Array<double, 1> &GetGrid()
	{
		BZPRECONDITION(Grid.extent(0) > 0);
		return Grid;
	}

	const blitz::Array<double, 1> &GetWeights()
	{
		BZPRECONDITION(Weights.extent(0) > 0);
		return Weights;
	}

	int GetGridIndex(const CoupledIndex &c)
	{
		BZPRECONDITION(IndexMap.count(c) == 1);
		CoupledIndexMap::iterator it = IndexMap.find(c);
		return it->second;
	}

	int IsGridIndex(const CoupledIndex &c)
	{
		BZPRECONDITION(IndexMap.count(c) < 2);
		CoupledIndexMap::iterator it = IndexMap.find(c);
		return (it != IndexMap.end()) ? it->second : -1;
	}

	CoupledIndex GetCoupledIndex(int i)
	{
		BZPRECONDITION(IndexList.extent(0) > i);
		return IndexList(i);
	}

	void BeginIndexList()
	{
		IndexMap.clear();
		curGridIndex = 0;
	}

	int AddIndex(const CoupledIndex &c)
	{
		IndexMap.insert( CoupledIndexMap::value_type(c, curGridIndex++) );
		return curGridIndex-1;
	}

	void EndIndexList()
	{
		//Resize arrays
		IndexList.resize(curGridIndex);
		Grid.resize(curGridIndex);
		Weights.resize(curGridIndex);

		typedef CoupledIndexMap::iterator Iterator;
		for (Iterator it=IndexMap.begin(); it!=IndexMap.end(); it++)
		{
			int i = (*it).second;
			BZPRECONDITION(i < curGridIndex);
		    IndexList(i) = (*it).first;
		}
		Grid = blitz::tensor::i;
		Weights = 1.0;

		curGridIndex = -1;
	}

};

} //Namespace

#endif


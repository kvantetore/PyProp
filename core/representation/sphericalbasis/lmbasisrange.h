#ifndef LMBASISRANGE_H
#define LMBASISRANGE_H

//#include <boost/unordered_map.hpp>
#include <map>
//#include <boost/hash.hpp>
//#include <boost/foreach.hpp>

#include "../../common.h"
#include "../orthogonalrepresentation.h"
#include <gsl/gsl_sf_coupling.h>


namespace SphericalBasis
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

class LmIndex
{
public:
	typedef int IndexType;

	IndexType l;
	IndexType m;

	LmIndex() : l(-1), m(-1) {}
	LmIndex(IndexType l, IndexType m) : l(l), m(m) {}

	bool operator==(LmIndex const& other) const
	{
		return l == other.l && m == other.m;
	}

	unsigned int GetInt() const
	{
		typedef unsigned int u;
		return ((u)l<<7) + m;
	}
	
};

class LmIndexCompareLT
{
public:
	bool operator()(LmIndex const& p1, LmIndex const& p2) const
	{
		return p1.GetInt() < p2.GetInt();
	}
};

//typedef boost::unordered_map<LmIndex, int> LmIndexMap;
typedef std::map<LmIndex, int, LmIndexCompareLT> LmIndexMap;

class LmBasisRange
{
private:
	blitz::Array<double, 1> Grid;		//Grid is a double representation of the index, it is only 
										//used to satisfy the Representation interface
	blitz::Array<double, 1> Weights;    //Weights is the Clebsch-Gordan coefficients 

	blitz::Array<LmIndex, 1> IndexList;  //Maps from a grid index to a lm index
	LmIndexMap IndexMap;				  //Maps from a lm index to a grid index

	int curGridIndex;
	
public:
	//Constructors
	LmBasisRange() : curGridIndex(0) {}
	~LmBasisRange() {}

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

	int GetGridIndex(const LmIndex &c)
	{
		BZPRECONDITION(IndexMap.count(c) == 1);
		LmIndexMap::iterator it = IndexMap.find(c);
		return it->second;
	}

	int IsGridIndex(const LmIndex &c)
	{
		BZPRECONDITION(IndexMap.count(c) < 2);
		LmIndexMap::iterator it = IndexMap.find(c);
		return (it != IndexMap.end()) ? it->second : -1;
	}

	LmIndex GetLmIndex(int i)
	{
		BZPRECONDITION(IndexList.extent(0) > i);
		return IndexList(i);
	}

	void BeginIndexList()
	{
		IndexMap.clear();
		curGridIndex = 0;
	}

	int AddIndex(const LmIndex &c)
	{
		IndexMap.insert( LmIndexMap::value_type(c, curGridIndex++) );
		return curGridIndex-1;
	}

	void EndIndexList()
	{
		//Resize arrays
		IndexList.resize(curGridIndex);
		Grid.resize(curGridIndex);
		Weights.resize(curGridIndex);

		typedef LmIndexMap::iterator Iterator;
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


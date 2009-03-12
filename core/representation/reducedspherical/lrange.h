#ifndef REDUCEDSPHERICAL_LRANGE_H
#define REDUCEDSPHERICAL_LRANGE_H

#include "../../common.h"

namespace ReducedSpherical
{

class LRange
{
private:
	blitz::Array<double, 1> Grid;
	blitz::Array<double, 1> Weights;
	
public:
	//Public fields
	int MaxL;
	int M;

	//Constructors
	LRange() : MaxL(0), M(0) {}
	LRange(int maxL) : MaxL(maxL), M(0) {}
	LRange(int maxL, int m) : MaxL(maxL), M(m) {}

	//Returns the number possible lm values
	inline int Count() const
	{
		return MaxL + 1 - std::abs(M);	
	}

	const blitz::Array<double, 1> &GetGrid()
	{
		if (Grid.size() == 0)
		{
			Grid.resize(Count());
			Grid = blitz::tensor::i;
			Grid += std::abs(M);
		}
		return Grid;
	}

	const blitz::Array<double, 1> &GetWeights()
	{
		if (Weights.size() == 0)
		{
			Weights.resize(Count());
			Weights = 1;
		}
		return Weights;
	}

};

} //Namespace

#endif

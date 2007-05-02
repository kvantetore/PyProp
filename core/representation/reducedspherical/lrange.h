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

	//Constructors
	LRange() : MaxL(0) {}
	LRange(int maxL) : MaxL(maxL) {}

	//Returns the number possible lm values
	inline int Count() const
	{
		return MaxL + 1;	
	}

	const blitz::Array<double, 1> &GetGrid()
	{
		if (Grid.size() == 0)
		{
			Grid.resize(Count());
			Grid = blitz::tensor::i;
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

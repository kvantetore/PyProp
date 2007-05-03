#ifndef ORTHOPOLRANGE_H
#define ORTHOPOLRANGE_H

#include "../../common.h"
#include "../../transform/orthopol/orthopoltools.h"

class OrthoPolRange
{
private:
	blitz::Array<double, 1> Grid;
	blitz::Array<double, 1> Weights;
	
public:
	int Count;
	OrthoPol::Parameter Param;
	double Alpha;
	
	//Constructors--------------------------------
	OrthoPolRange() :
		Count(0),
		Alpha(0.0)
	{}

	OrthoPolRange(const OrthoPol::Parameter &param, int N) :
		Count(N),
		Param(param),
		Alpha(0.0)
	{}
	
	//Member functions----------------------------
	blitz::Array<double, 1>& GetGrid()
	{
		if (Grid.extent(0) == 0)
		{
			OrthoPol::GridAndWeights(Count, Param, Grid, Weights, Alpha);
		}
		return Grid;
	}

	blitz::Array<double, 1>& GetWeights()
	{
		if (Weights.extent(0) == 0)
		{
			OrthoPol::GridAndWeights(Count, Param, Grid, Weights, Alpha);
		}
		return Weights;
	}	
	
};


#endif


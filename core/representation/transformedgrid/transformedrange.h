#ifndef TRANSFORMEDRANGE_H
#define TRANSFORMEDRANGE_H

#include "../../common.h"
#include "../../transform/transformedgrid/tools.h"

class TransformedRange
{
private:
	blitz::Array<double, 1> Grid;
	blitz::Array<double, 1> Weights;
	
public:
	int Count;
	TransformedGrid::Parameter Param;
	
	//Constructors--------------------------------
	TransformedRange() :
		Count(0)
	{}

	TransformedRange(const TransformedGrid::Parameter &param, int N) :
		Count(N),
		Param(param)
	{}
	
	//Member functions----------------------------
	blitz::Array<double, 1>& GetGrid()
	{
		if (Grid.extent(0) == 0)
		{
			blitz::Array<double, 1> X, Y;
			TransformedGrid::XYmat(Count, Param, X, Y, Grid);
		}
		return Grid;
	}

	blitz::Array<double, 1>& GetWeights()
	{
		if (Weights.extent(0) == 0)
		{
			TransformedGrid::SetupWeights(Count, Param, GetGrid(), Weights);
		}
		return Weights;
	}
	
};


#endif


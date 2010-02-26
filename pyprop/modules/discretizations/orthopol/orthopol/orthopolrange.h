#ifndef ORTHOPOLRANGE_H
#define ORTHOPOLRANGE_H

#include "../../common.h"
#include "../../transform/orthopol/orthopoltools.h"

namespace OrthoPol
{

class OrthoPolRange
{
private:
	blitz::Array<double, 1> Grid;
	blitz::Array<double, 1> Weights;
	
public:
	int Count;
	OrthoPol::Parameter Param;
	
	//Constructors--------------------------------
	OrthoPolRange() :
		Count(0)
	{}

	OrthoPolRange(const Parameter &param, int N) 
	{
		Initialize(param, N);
	}
	
	//Member functions----------------------------
	blitz::Array<double, 1>& GetGrid()
	{
		return Grid;
	}

	blitz::Array<double, 1>& GetWeights()
	{
		return Weights;
	}

	void Initialize(const Parameter &param, int n)
	{
		ScaledGridAndWeights(n, param, Grid, Weights);
		Count = n;
		Param = param;
	}
	
};

} //Namespace

#endif


#ifndef REDUCEDSPHERICAL_OMEGARANGE_H
#define REDUCEDSPHERICAL_OMEGARANGE_H

#include "../../common.h"
#include "../../transform/reducedspherical/reducedsphericaltools.h"

namespace ReducedSpherical
{

class ThetaRange
{
private:
	//This grid represents the point (theta) at the given index, and
	//is what really should be used in most cases
	blitz::Array<double, 1> Grid;

	// Array of the weights used in theta integration
	blitz::Array<double, 1> Weights;

public:
	int MaxL;

	//Constructors
	ThetaRange() {}

	void SetupRange(int maxL)
	{
		MaxL = maxL;
		
		ReducedSphericalTools trans;
		trans.Initialize(maxL);

		Grid.resize(trans.GetThetaGrid().shape());
		Grid = trans.GetThetaGrid();

		Weights.resize(trans.GetWeights().shape());
		Weights = trans.GetWeights();
	}

	int Count()
	{
		return Grid.size();
	}

	const blitz::Array<double, 1> &GetGrid()
	{
		return Grid;
	}

	const blitz::Array<double, 1> &GetWeights()
	{
		return Weights;
	}	
};

} //namespace

#endif


#ifndef OMEGARANGE_H
#define OMEGARANGE_H

#include "../../common.h"
#include "../../transform/spherical/shtools.h"
#include <iostream>

class OmegaRange
{
private:
	//This array is just to make it easier to comply with the 
	//Representation interface, it should probably not be used.
	//It contains only the index value, stored as a double!
	blitz::Array<double, 1> IndexGrid;

	//This grid represents the point (theta, phi) at the given index, and
	//is what really should be used in most cases
	//OmegaGrid(:, 0) = theta(:)
	//OmegaGrid(:, 1) = phi(:)
	blitz::Array<double, 2> OmegaGrid;

	// Array of the weights used in theta integration
	blitz::Array<double, 1> Weights;

public:
	enum GridType
	{
		SloanWomersley,
		Equidistant,
		Gauss
	};

	int Type;
	int MaxL;

	//Constructors
	OmegaRange() {}

	void SetupRange(GridType type, int maxL)
	{
		Type = type;
		MaxL = maxL;

		if (type != Gauss)
		{
			throw std::runtime_error("Invalid Spherical Grid Type. Only Gauss is currently supported");
		}

		SphericalTransformTensorGrid trans;
		trans.Initialize(maxL);

		OmegaGrid.resize(trans.GetOmegaGrid().shape());
		OmegaGrid = trans.GetOmegaGrid();

		Weights.resize(trans.GetWeights().shape());
		Weights = trans.GetWeights();
	}

	int Count()
	{
		return OmegaGrid.size();
	}

	const blitz::Array<double, 1> &GetIndexGrid()
	{
		if (IndexGrid.size() == 0)
		{
			std::cout << "Warning: Calling GetIndexGrid() on OmegaRange for the first time. "
			          << "You should probably not use this grid, as it is only a list of "
				  << "the index values which is propably not what you are expecting. "
				  << "Consider writing a specialized loop instead. See SphericalDynamicPotentialEvaluator "
				  << std::endl;
			IndexGrid.resize(Count());
			IndexGrid = blitz::tensor::i;
		}
		return IndexGrid;
	}	
	
	const blitz::Array<double, 2> &GetOmegaGrid()
	{
		return OmegaGrid;
	}

	const blitz::Array<double, 1> &GetWeights()
	{
		return Weights;
	}	
};

#endif


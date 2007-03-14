#ifndef CARTESIANRANGE_H
#define CARTESIANRANGE_H

#include "../common.h"

class CartesianRange
{
private:
	blitz::Array<double, 1> Grid;
	blitz::Array<double, 1> Weights;
	
public:
	double Min;
	double Max;
	double Dx;
	int Count;
	
	//If TranslatedGrid == true, this class will create a grid
	// [0:Max Min:0]
	bool TranslatedGrid;
	
	//Constructors--------------------------------
	CartesianRange() :
		Min(0),
		Max(1),
		Dx(1),
		Count(1)
	{}
	
	CartesianRange(double min, double max, int count, bool translatedGrid=false) :
		Min(min),
		Max(max),
		Count(count),
		TranslatedGrid(translatedGrid)
	{
		Dx = (max - min) / (double)count;
	}
	
	CartesianRange(const CartesianRange& r2)
	{
		*this = r2;
	}
	
	//Member functions----------------------------
	blitz::Array<double, 1>& GetGrid()
	{
		if (Grid.extent(0) == 0)
		{
			Grid.resize(Count);
			if (TranslatedGrid)
			{
				Grid = Min + Dx * ((blitz::tensor::i + (Count/2)) % Count) ;
			} 
			else
			{
				Grid = Min + Dx * blitz::tensor::i ;
			}
			
		}
		return Grid;
	}

	blitz::Array<double, 1>& GetWeights()
	{
		if (Weights.extent(0) == 0)
		{
			Weights.resize(Count);
			Weights = Dx;
		}
		return Weights;
	}
	
	double GetPosition(int index)
	{
		return Min + index * Dx;
	}
	
	int GetIndex(double position)
	{
		return (int)round((position - Min) / Dx);
	}
	
	CartesianRange& operator=(const CartesianRange &r2)
	{
		Min = r2.Min;
		Max = r2.Max;
		Count = r2.Count;
		Dx = r2.Dx;
		TranslatedGrid = r2.TranslatedGrid;
		Grid.resize(0);
		
		return *this;
	}
};


#endif

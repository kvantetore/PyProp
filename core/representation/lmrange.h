#ifndef LMRANGE_H
#define LMRANGE_H

#include "../common.h"

class LmRange
{
private:
	//This array is just to make it easier to comply with the 
	//Representation interface, it should probably not be used.
	//It contains only the index value, stored as a double!
	blitz::Array<double, 1> IndexGrid;
	blitz::Array<double, 2> LmGrid;
	blitz::Array<double, 1> Weights;
	
public:
	//Public fields
	int MaxL;

	//Constructors
	LmRange() : MaxL(0) {}
	LmRange(int maxL) : MaxL(maxL) {}

	//Returns the number possible lm values
	inline int Count() const
	{
		return sqr(MaxL + 1);	
	}

	const blitz::Array<double, 1> &GetIndexGrid()
	{
		if (IndexGrid.size() == 0)
		{
			std::cout << "Warning: Calling GetIndexGrid() on LmRange for the first time. "
			          << "You should probably not use this grid, as it is only a list of "
				  << "the index values which is propably not what you are expecting. "
				  << "Consider writing a specialized loop instead "
				  << std::endl;
			IndexGrid.resize(Count());
			IndexGrid = blitz::tensor::i;
		}
		return IndexGrid;
	}

	const blitz::Array<double, 2> &GetLmGrid()
	{
		if (LmGrid.size() == 0)
		{
			LmGrid.resize(Count(),2);
			for (int i=0; i<Count(); i++)
			{
				LmGrid(i, 0) = GetL(i);
				LmGrid(i, 1) = GetM(i);
			}
		}
		return LmGrid;
	}
	
	blitz::Array<double, 1> GetWeights()
	{
		if (Weights.size() == 0)
		{
			Weights.resize(Count());
		}
		return Weights;
	}

	/** -------------- Mapping between index and l,m pairs ----------------- */

	//Maps (l,m) to the linear index
	static inline int MapLmIndex(int l, int m)
	{
		return (m + l) + l*l;
	}

	//Gets the l-value from the linear index
	static inline int GetL(int index)
	{
		return static_cast<int>(sqrt((double)index));		
	}

	//Gets the m-value from the linear index
	static inline int GetM(int index)
	{
		int l = GetL(index);
		return index - l - sqr(l);
	}
	
	//Gets both l and m value from the linear index as a tinyvector of length 2
	//the first index is the l value, and the second index is the m value
	static inline blitz::TinyVector<int, 2> MapLmIndexInverse(int index)
	{
		blitz::TinyVector<int,2> lm;
		int l = GetL(index);
		lm(0) = l;
		lm(1) = index - l - sqr(l);
		return lm;
	}

};

#endif

#ifndef BLITZUTILS_H
#define BLITZUTILS_H

#include "../common.h"

template<int Rank>
int inline CoordToIndex(const blitz::TinyVector<int, Rank> &shape, const blitz::TinyVector<int, Rank> &coord)
{
	int index = coord(0);
	for (int i=1; i<Rank; i++)
	{
		index = (index*shape(i)) + coord(i);	
	}
	return index;
}

int inline CoordToIndex(const blitz::Array<int, 1> &shape, const blitz::Array<int, 1> &coord)
{
	int index = coord(0);
	for (int i=1; i<shape.size(); i++)
	{
		index = (index*shape(i)) + coord(i);
	}
	return index;
}

template<int Rank>
blitz::TinyVector<int, Rank> inline IndexToCoord(const blitz::TinyVector<int, Rank> &shape, int index)
{
	blitz::TinyVector<int, Rank> coord;
	for (int i=0; i<Rank; i++)
	{
		coord(i) = index % shape(i);
		index -= coord(i) * shape(i);
	}
	return coord;
}

blitz::Array<int, 1> inline IndexToCoord(const blitz::Array<int, 1> &shape, int index)
{
	blitz::Array<int, 1> coord(shape.size());
	for (int i=0; i<shape.size(); i++)
	{
		coord(i) = index % shape(i);
		index -= coord(i) * shape(i);
	}
	return coord;
}

template<class T, int Rank>
bool Contains(const blitz::Array<T, Rank> &vector, int value)
{
	for (int i=0; i<vector.size(); i++)
	{
		if (vector(i) == value)
		{
			return true;
		}
	}
	return false;
}

#endif


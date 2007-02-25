#ifndef INDEX_H
#define INDEX_H

#include "common.h"

template<int Rank>
inline bool IncIndex(const blitz::TinyVector<int, Rank>& extent, blitz::TinyVector<int, Rank> &index, int curDimension, int skipDimension)
{
	if (curDimension < 0) {
		return false;
	}
	
	if (skipDimension == curDimension)
	{
		return IncIndex(extent, index, curDimension - 1, skipDimension);
	}
	
	index(curDimension)++;
	if (index(curDimension) >= extent(curDimension))
	{
		index(curDimension) = 0;
		return IncIndex(extent, index, curDimension - 1, skipDimension);
	}
	
	return true;
}

#endif


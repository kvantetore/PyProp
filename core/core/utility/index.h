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

/* Example usage
 
   typedef blitz::TinyVector<int, Rank> IndexType

   IndexType startIndex = 0;     //The current index
   IndexTyp  endIndex = 0;
   int skipDimension = r;   //The dimension (rank) which we should not iterate over

   

   do 
   {
   	   endIndex = startIndex;
	   endIndex(skipDimension) = data.extent(skipDimension) - 1;
       blitz::Range
       blitz::Array<cplx, N> = slice;
   } while (IncIndex(data.shape(), index, N-1, skipDimension)


*/

#endif


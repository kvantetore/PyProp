#ifndef COMPRESSEDREPRESENTATION_H
#define COMPRESSEDREPRESENTATION_H

#include "orthogonalrepresentation.h"

/*
 * CompressedRepresentation is a representation, where the data storage is
 * one dimensional, but represents a higher dimensional physical grid.
 *
 * An example of this, is the angular and spherical harmonic representation 
 * in spherical coordinates, where theta, phi and l, m is compressed to 
 * one rank, but the potentials need the expaned grid values.
 *
 * Get...ExpandedGrid() returns a 2D array, where the first rank is
 * the grid point index, and the second rank is the expanded dimension index
 */
class CompressedRepresentation : public OrthogonalRepresentation
{
public:
	typedef shared_ptr< CompressedRepresentation > Ptr;
	typedef blitz::Array<double, 2> GridArray;

	virtual GridArray GetGlobalExpandedGrid() = 0;

	GridArray GetLocalExpandedGrid()
	{
		GridArray globalGrid( this->GetGlobalExpandedGrid() );
		int size = globalGrid.extent(0);

		blitz::Range indexRange = this->GetDistributedModel()->GetLocalIndexRange(size, GetBaseRank());
		return globalGrid(indexRange, blitz::Range::all());
	}
};


#endif


#ifndef CUSTOMGRIDREPRESENTATION_H
#define CUSTOMGRIDREPRESENTATION_H

#include "../common.h"
#include "../mpi/distributedmodel.h"
#include "representation.h"
#include "orthogonalrepresentation.h"

class CustomGridRepresentation : public OrthogonalRepresentation
{
public:
	typedef boost::shared_ptr< CustomGridRepresentation > Ptr;
	
	/*
	 * Constructors
	 */
	CustomGridRepresentation() {}
	
	virtual Representation<1>::RepresentationPtr Copy()
	{
		return Representation<1>::RepresentationPtr(new CustomGridRepresentation(*this));
	}
	
	/* 
	 * Returns the global (distributed) grid 
	 */
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank)
	{
		int effectiveRank = rank - this->GetBaseRank();
		if (effectiveRank != 0)
		{
			cout << "WARNING: specified invalid rank " << rank << " for 1d rank " << this->GetBaseRank() << endl;
		}
		return GlobalGrid;
	}	

	/* 
	 * Returns the global (distributed) weights
	*/
	virtual blitz::Array<double, 1> GetGlobalWeights(int rank)
	{
		int effectiveRank = rank - this->GetBaseRank();
		if (effectiveRank != 0)
		{
			cout << "WARNING: specified invalid rank " << rank << " for 1d rank " << this->GetBaseRank() << endl;
		}
		return GlobalWeights;
	}	

	/*
	 * Returns the number of grid points
	 */
	virtual blitz::TinyVector<int, 1> GetFullShape()
	{
		return blitz::TinyVector<int, 1>(GlobalGrid.size());
	}

	virtual void ApplyConfigSection(const ConfigSection &cfg);

private:
	blitz::Array<double, 1> GlobalGrid;
	blitz::Array<double, 1> GlobalWeights;
};

#endif

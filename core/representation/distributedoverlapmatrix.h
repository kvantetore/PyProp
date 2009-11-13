#ifndef DISTRIBUTEDOVERLAPMATRIX_H
#define DISTRIBUTEDOVERLAPMATRIX_H

#include "../wavefunction.h"
#include <boost/shared_ptr.hpp>

template<int Rank>
class DistributedOverlapMatrix
{
public:
	DistributedOverlapMatrix() {}
	virtual ~DistributedOverlapMatrix() {}
	
	virtual void MultiplyOverlapRank(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &destPsi, int opRank);
};

#endif


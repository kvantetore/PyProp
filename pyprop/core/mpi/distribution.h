#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "../common.h"

/* 
 * Class to keep track of the current distribution among processors.
 */
class Distribution
{
public:
	typedef shared_ptr< Distribution > Ptr;
	typedef blitz::Array<int, 1> DataArray;

	Distribution() : Distrib(0) {}

	Distribution(int procRank)
	{
		Distrib.resize(procRank);
	}

	Distribution(const Distribution &other)
	{
		if (other.Distrib.size() > 0)
		{
			this->Distrib.resize(other.Distrib.shape());
			this->Distrib = other.Distrib;
		}
	}

	int GetProcRank() const
	{
		return Distrib.extent(0);
	}
		
	const DataArray GetDistribution() const
	{
		return Distrib;
	}

	void SetDistribution(const DataArray &distrib)
	{
		Distrib = distrib;
	}

private:
	DataArray Distrib;
};
typedef boost::shared_ptr<Distribution> DistributionPtr;


#endif


#ifndef COMBINEDABSORBER_H
#define COMBINEDABSORBER_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/representation.h"
#include "../utility/blitztricks.h"

/*
 * Absorbing potential of the form cos(r)**80
 */

class AbsorberModel
{
public:
	typedef shared_ptr<AbsorberModel> Ptr;
	typedef blitz::Array<double, 1> Vector;

	virtual ~AbsorberModel() {}
	virtual void ApplyConfigSection(const ConfigSection &config) {}
	virtual void SetupStep(Vector grid) = 0;
	virtual Vector GetScaling() = 0;
};

template<int Rank>
class CombinedAbsorberPotential
{
private:
	typedef std::pair< AbsorberModel::Ptr, int> AbsorberRankElement;
	std::vector< AbsorberRankElement > AbsorberList;
		
public:
	void AddAbsorber(AbsorberModel::Ptr absorber, int rank)
	{
		AbsorberList.push_back( AbsorberRankElement(absorber, rank) );
	}
	
	/** Updates the wavefunction a step with this dynamic potential **/
	void ApplyPotential(Wavefunction<Rank> &psi, cplx timeStep, double curTime)
	{
		//Set up grid
		typename Representation<Rank>::Ptr repr = psi.GetRepresentation();

		//Iterate
		typedef typename blitz::Array<cplx, Rank>::iterator iterator;

		for (unsigned int absorberIndex=0; absorberIndex<AbsorberList.size(); absorberIndex++)
		{
			blitz::Array<double, 1> scaling = AbsorberList[absorberIndex].first->GetScaling();
			int rank = AbsorberList[absorberIndex].second;

			blitz::Array<cplx, 3> data = MapToRank3(psi.Data, rank, 1);
			cplx* ptr = data.data();
			for (int outerIndex=0; outerIndex<data.extent(0); outerIndex++)
			{
				for (int index=0; index<data.extent(1); index++)
				{
					double scalingValue = scaling(index);
					for(int innerIndex=0; innerIndex<data.extent(2); innerIndex++)
					{
						*ptr *= scalingValue;
						ptr++;
					}
				}
			}
			
			/*
			for (iterator it=psi.Data.begin(); it!=psi.Data.end(); it++)
			{
				(*it) *= scaling(it.position()(rank));
			}
			*/
		}
	}
};


#endif


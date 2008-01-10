#ifndef SPHERICALABSORBER_H
#define SPHERICALABSORBER_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/sphericalrepresentation.h"
#include "../representation/radialrepresentation.h"

/*
 * Absorbing potential of the form cos(r)**80
 */

class AbsorberModel
{
private:
	blitz::Array<double, 1> Scaling;
		
public:
	void ApplyConfigSection(const ConfigSection &config, int rank)
	{
		//TODO: Support different absorber models
	}

	const blitz::Array<double, 1>& GetScaling(const blitz::Array<double, 1> &grid, double gridmax)
	{
		if (Scaling.size() != grid.size())
		{
			Scaling.resize(grid.size());
			Scaling = 1.0 - pow(cos(M_PI/2.0 * (1.0 - fabs(grid) / gridmax)), 80);
		}

		return Scaling;
	}
};

template<int Rank>
class SphericalAbsorbingPotential
{
private:
	blitz::TinyVector<AbsorberModel, Rank-1> Absorber;
		
public:
	
	void ApplyConfigSection(const ConfigSection &config)
	{
		for (int i=0; i<Rank-1; i++)
		{
			Absorber(i).ApplyConfigSection(config, i);
		}
	}

	/** Updates the wavefunction a step with this dynamic potential **/
	void ApplyPotential(Wavefunction<Rank> &psi, cplx timeStep, double curTime)
	{
		//Set up grid
		typename CombinedRepresentation<Rank>::Ptr repr = dynamic_pointer_cast< CombinedRepresentation<Rank> >(psi.GetRepresentation());

		blitz::TinyVector< blitz::Array<double, 1>, Rank-1 > grid;
		blitz::TinyVector< blitz::Array<double, 1>, Rank-1 > scaling;
		for (int i=0; i<Rank-1; i++)
		{
			//grid
			grid(i).reference(repr->GetLocalGrid(i));
			blitz::Array<double, 1> globalGrid(repr->GetGlobalGrid(i));

			//Max value
			//TODO: Do properly
			int gridSize = globalGrid.size();
			double maxR = max(globalGrid(gridSize-1), globalGrid(0));
			
			//Setup absorber
			scaling(i).reference( Absorber(i).GetScaling(grid(i), maxR) );
		}

		//Iterate
		typename blitz::Array<cplx, Rank>::iterator it = psi.Data.begin();
		double curScale = 1;
		for (int linearCount=0; linearCount<psi.Data.size(); linearCount++)
		{
			curScale = 1;
			for (int i=0; i<Rank-1; i++)
			{
				curScale *= scaling(i)(it.position()(i));
			}
			
			(*it) *= curScale;
			
			it++;
		}
	}
};


#endif


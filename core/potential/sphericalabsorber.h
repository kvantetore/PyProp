#ifndef SPHERICALABSORBER_H
#define SPHERICALABSORBER_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/sphericalrepresentation3d.h"
#include "../representation/radialrepresentation.h"

/*
 * Absorbing potential of the form cos(r)**80
 */
template<int Rank>
class SphericalAbsorbingPotential
{
public:
	/** Updates propagates the wavefunction a step with this dynamic potential **/
	void ApplyPotential(Wavefunction<Rank> &psi, cplx timeStep, double curTime)
	{
		//Set up grid
		blitz::Array<double, 1> radialGrid;
		SphericalRepresentation3D *repr = static_cast<SphericalRepresentation3D*>(&psi.GetRepresentation());
		radialGrid.reference(repr->GetLocalGrid(0));

		int radialCount = psi.Data.extent(0);
		int omegaCount = psi.Data.extent(1);
		if (radialCount * omegaCount != psi.Data.size())
		{
			throw std::runtime_error("Invalid wavefunction size in SphericalAbsorbingPotential");
		}

		//Max value
		double maxR = radialGrid(radialGrid.extent(0)-1);
		
		//Iterate
		blitz::TinyVector<double, Rank> pos;
		blitz::TinyVector<int, Rank> indexPos;
		typename blitz::Array<cplx, Rank>::iterator it = psi.Data.begin();

		//WARNING: this iteration assumes that the radial dim has max stride!
		for (int rIdx=0; rIdx<radialCount; rIdx++)
		{
			double r = radialGrid(rIdx);
			double scale = 1.0 - pow(cos(M_PI/2.0 * (1.0 - fabs(r) / maxR)), 80);
			
			for (int omegaIdx=0; omegaIdx<omegaCount; omegaIdx++)
			{
				*it *= scale;
				it++;
			}
		}
	}
};


#endif


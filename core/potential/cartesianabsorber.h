#ifndef CARTESIANABSORBER_H
#define CARTESIANABSORBER_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/cartesianrepresentation.h"

/*
 * Absorbing potential of the form cos()**80
 */
template<int Rank>
class CartesianAbsorbingPotential
{
public:
	/** Updates propagates the wavefunction a step with this dynamic potential **/
	void ApplyPotential(Wavefunction<Rank> &psi, cplx timeStep, double curTime)
	{
		//Set up grid
		blitz::TinyVector< blitz::Array<double, 1>, Rank> grid;
		Representation<Rank> *reprBasic = &psi.GetRepresentation();
		CartesianRepresentation<Rank> *repr = (CartesianRepresentation<Rank>*) reprBasic;
		for (int curRank = 0; curRank<Rank; curRank++)
		{
			grid(curRank).reference(repr->GetLocalGrid(curRank));
		}

		//Max values
		blitz::TinyVector<double, Rank> maxPos;
		for (int curRank = 0; curRank<Rank; curRank++)
		{
			maxPos(curRank) = repr->GetRange(curRank).Max;
		}
		
		//Iterate
		blitz::TinyVector<double, Rank> pos;
		blitz::TinyVector<int, Rank> indexPos;
		typename blitz::Array<cplx, Rank>::iterator it = psi.Data.begin();
		for (int linearCount=0; linearCount<psi.Data.size(); linearCount++)
		{
			for (int curRank=0; curRank<Rank; curRank++)
			{
				pos(curRank) = grid(curRank)(it.position()(curRank));
			}
			
			double scale = 1;
			for (int curRank=0; curRank<Rank; curRank++)
			{ 
				scale *= 1.0 - pow(cos(M_PI/2.0 * (1.0 - fabs(pos(curRank)) / maxPos(curRank))), 80);
			}
			*it *= scale;
			
			it++;
		}
	}
};


#endif


#ifndef WAVEFUNCTIONEVALUATOR_H
#define WAVEFUNCTIONEVALUATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/representation.h"

template<class InitialConditionClass, int Rank>
class WavefunctionEvaluator : public InitialConditionClass
{
public:
	/** Setup the wavefunction **/
	void SetupWavefunction(Wavefunction<Rank> &psi)
	{
		//Set up grid
		blitz::TinyVector< blitz::Array<double, 1>, Rank> grid;
		typename Representation<Rank>::Ptr repr = psi.GetRepresentation();
		for (int curRank = 0; curRank<Rank; curRank++)
		{
			grid(curRank).reference(repr->GetLocalGrid(curRank));
		}
	
		//Iterate through each point in the wavefunction
		blitz::TinyVector<double, Rank> pos;
		blitz::TinyVector<int, Rank> indexPos;
		typename blitz::Array<cplx, Rank>::iterator it = psi.Data.begin();
		for (int linearCount=0; linearCount<psi.Data.size(); linearCount++)
		{
			for (int curRank=0; curRank<Rank; curRank++)
			{
				pos(curRank) = grid(curRank)(it.position()(curRank));
			}
			
			*it = GetInitialValue(pos);
			
			it++;
		}

		psi.Normalize();
	}

};

#endif


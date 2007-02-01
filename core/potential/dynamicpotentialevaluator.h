#ifndef DYNAMICPOTENTIALEVALUATOR_H
#define DYNAMICPOTENTIALEVALUATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/representation.h"
#include "staticpotential.h"
#include "examplepotentials.h"

template<class DynamicPotentialClass, int Rank>
class DynamicPotentialEvaluator : public DynamicPotentialClass
{
public:
	/** Updates propagates the wavefunction a step with this dynamic potential **/
	void ApplyPotential(Wavefunction<Rank> &psi, cplx timeStep, double curTime)
	{
		//Set up PotentialClass
		this->CurTime = curTime;
		this->TimeStep = timeStep;
	
		//Set up grid
		blitz::TinyVector< blitz::Array<double, 1>, Rank> grid;
		Representation<Rank> *repr = &psi.GetRepresentation();
		for (int curRank = 0; curRank<Rank; curRank++)
		{
			grid(curRank).reference(repr->GetLocalGrid(psi, curRank));
		}
		
		blitz::TinyVector<double, Rank> pos;
		blitz::TinyVector<int, Rank> indexPos;
		typename blitz::Array<cplx, Rank>::iterator it = psi.Data.begin();
		for (int linearCount=0; linearCount<psi.Data.size(); linearCount++)
		{
			for (int curRank=0; curRank<Rank; curRank++)
			{
				pos(curRank) = grid(curRank)(it.position()(curRank));
			}
			
			*it *= exp( - I * GetPotentialValue(pos) * timeStep );
			
			it++;
		}
	}

	/** Updates a static potential with the expotential of the potential of this dynamic potential. */
	void UpdateStaticPotential(StaticPotential<Rank> &potential, const Wavefunction<Rank> &psi, cplx timeStep, double curTime)
	{
		//Set up PotentialClass
		this->CurTime = curTime;
		this->TimeStep = timeStep;
	
		//Set up grid
		blitz::TinyVector< blitz::Array<double, 1>, Rank> grid;
		Representation<Rank> *repr = &psi.GetRepresentation();
		for (int curRank = 0; curRank<Rank; curRank++)
		{
			grid(curRank).reference(repr->GetLocalGrid(psi, curRank));
		}
		
		blitz::TinyVector<double, Rank> pos;
		blitz::TinyVector<int, Rank> indexPos;
		typename blitz::Array<cplx, Rank>::iterator it = potential.GetPotentialData().begin();
		for (int linearCount=0; linearCount<psi.Data.size(); linearCount++)
		{
			for (int curRank=0; curRank<Rank; curRank++)
			{
				pos(curRank) = grid(curRank)(it.position()(curRank));
			}
			
			*it = exp( - I * GetPotentialValue(pos) * timeStep );
			
			it++;
		}	
	}

	/** 
   	  * Returns an array containing this dynamic potential at the given time. 
	  * WARNING: This method allocates and returns a new blitz array in the same size as the wavefunction
	  */
	blitz::Array<cplx, Rank> GetPotential(const Wavefunction<Rank> &psi, cplx timeStep, double curTime)
	{
		//Set up PotentialClass
		this->CurTime = curTime;
		this->TimeStep = timeStep;
	
		//Set up grid
		blitz::TinyVector< blitz::Array<double, 1>, Rank> grid;
		Representation<Rank> *repr = &psi.GetRepresentation();
		for (int curRank = 0; curRank<Rank; curRank++)
		{
			grid(curRank).reference(repr->GetLocalGrid(psi, curRank));
		}
		
		blitz::TinyVector<double, Rank> pos;
		blitz::TinyVector<int, Rank> indexPos;
		
		blitz::Array<cplx, Rank> potentialData(psi.Data.shape());

		typename blitz::Array<cplx, Rank>::iterator it = potentialData.begin();
		for (int linearCount=0; linearCount<potentialData.size(); linearCount++)
		{
			for (int curRank=0; curRank<Rank; curRank++)
			{
				pos(curRank) = grid(curRank)(it.position()(curRank));
			}
			
			*it = GetPotentialValue(pos);
			it++;
		}

		return potentialData;
	}

};

#endif


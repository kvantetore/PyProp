#ifndef DYNAMICPOTENTIALEVALUATOR_H
#define DYNAMICPOTENTIALEVALUATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/representation.h"
#include "staticpotential.h"
#include "potentialaction.h"
#include "potentialbase.h"


/* The class inherits from DynamicPotentialClass and ActionClass
 *   - DynamicPotentialClass defines the potential
 *   - ActionClass defines the how the potential is to be applied to the data
 */
template<class DynamicPotentialClass, class ActionClass, int Rank>
class DynamicPotentialEvaluatorBase : public DynamicPotentialClass, public ActionClass
{
public:
	/** Updates propagates the wavefunction a step with this dynamic potential **/
	void IterateAction(const Wavefunction<Rank> &psi, blitz::Array<cplx, Rank> updateData, const cplx &timeStep, const double &curTime)
	{
		//Set up PotentialClass
		this->CurTime = curTime;
		this->TimeStep = timeStep;

		//Get representations
		Representation<Rank> *repr = &psi.GetRepresentation();

		blitz::TinyVector< blitz::Array<double, 1>, Rank > grid;
		for (int i=0; i<Rank; i++)
		{
			grid(i).reference(repr->GetLocalGrid(i));
		}
		
		blitz::TinyVector<double, Rank> pos;
		typename blitz::Array<cplx, Rank>::iterator it = updateData.begin();
		for (int linearCount=0; linearCount<updateData.size(); linearCount++)
		{
			for (int i=0; i<Rank; i++)
			{
				pos(i) = grid(i)(it.position()(i));
			}
			
			// Uses the function from the inherited classes
			ApplyAction(it, GetPotentialValue(pos), timeStep );
			it++;
		}
	}
};

template <class PotentialClass, int Rank>
class DynamicPotentialEvaluator
{
public:
	typedef ApplyPotentialClass<Rank> ApplyClass;
	typedef UpdatePotentialClass<Rank> UpdateClass;
	typedef GetPotentialClass<Rank> GetClass;
		
	DynamicPotentialEvaluatorBase< PotentialClass, ApplyClass, Rank> Apply;
	DynamicPotentialEvaluatorBase< PotentialClass, UpdateClass, Rank> Update;
	DynamicPotentialEvaluatorBase< PotentialClass, GetClass, Rank> Get;

	
	void ApplyConfigSection(const ConfigSection &config)
	{
		Apply.ApplyConfigSection(config);
		Update.ApplyConfigSection(config);
		Get.ApplyConfigSection(config);
	}

	/** Updates propagates the wavefunction a step with this dynamic potential **/
	void ApplyPotential(Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		Apply.IterateAction(psi, psi.GetData(), timeStep, curTime);
	}

	/** Updates a static potential with the expotential of the potential of this dynamic potential. */
	void UpdateStaticPotential(StaticPotential<Rank> &potential, const Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		blitz::Array<cplx, Rank> potentialData( potential.GetPotentialData() );
		Update.IterateAction(psi, potentialData, timeStep, curTime);
	}

	/** 
   	  * Returns an array containing this dynamic potential at the given time. 
	  * WARNING: This method allocates and returns a new blitz array in the same size as the wavefunction
	  */
	blitz::Array<cplx, Rank> GetPotential(const Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		//TODO: Fix for multiproc
		blitz::Array<cplx, Rank> potentialData(psi.Data.shape());
		Get.IterateAction(psi, potentialData, timeStep, curTime);
		return potentialData;
	}
};

#endif


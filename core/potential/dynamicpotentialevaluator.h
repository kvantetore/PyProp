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
template<class ActionClass, int Rank>
class DynamicPotentialEvaluatorBase : public ActionClass
{
public:
	/** Updates propagates the wavefunction a step with this dynamic potential **/
	template<class DynamicPotentialClass>
	void IterateAction(DynamicPotentialClass &potential, const Wavefunction<Rank> &psi, blitz::Array<cplx, Rank> updateData, const cplx &timeStep, const double &curTime)
	{
		//Set up PotentialClass
		potential.CurTime = curTime;
		potential.TimeStep = timeStep;

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
			ApplyAction(it, potential.GetPotentialValue(pos), timeStep );
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
	
	PotentialClass Potential;
	DynamicPotentialEvaluatorBase< ApplyClass, Rank> Apply;
	DynamicPotentialEvaluatorBase< UpdateClass, Rank> Update;
	DynamicPotentialEvaluatorBase< GetClass, Rank> Get;
	
	void ApplyConfigSection(const ConfigSection &config)
	{
		Potential.ApplyConfigSection(config);
	}

	/** Updates propagates the wavefunction a step with this dynamic potential **/
	void ApplyPotential(Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		Apply.IterateAction(Potential, psi, psi.GetData(), timeStep, curTime);
	}

	/** Updates a static potential with the expotential of the potential of this dynamic potential. */
	void UpdateStaticPotential(StaticPotential<Rank> &potential, const Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		blitz::Array<cplx, Rank> potentialData( potential.GetPotentialData() );
		Update.IterateAction(Potential, psi, potentialData, timeStep, curTime);
	}

	/** 
   	  * Returns an array containing this dynamic potential at the given time. 
	  * WARNING: This method allocates and returns a new blitz array in the same size as the wavefunction
	  */
	blitz::Array<cplx, Rank> GetPotential(const Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		//TODO: Fix for multiproc
		blitz::Array<cplx, Rank> potentialData(psi.Data.shape());
		Get.IterateAction(Potential, psi, potentialData, timeStep, curTime);
		return potentialData;
	}

	double CalculateExpectationValue(const Wavefunction<Rank> &psi, const cplx &timeStep, const double curTime)
	{
		blitz::Array<cplx, Rank> data(psi.Data);
			
		//Set up PotentialClass
		Potential.CurTime = curTime;
		Potential.TimeStep = timeStep;

		//Get representations
		Representation<Rank> *repr = &psi.GetRepresentation();

		blitz::TinyVector< blitz::Array<double, 1>, Rank > grid;
		for (int i=0; i<Rank; i++)
		{
			grid(i).reference(repr->GetLocalGrid(i));
		}

		blitz::TinyVector< blitz::Array<double, 1>, Rank > weights;
		for (int i=0; i<Rank; i++)
		{
			weights(i).reference(repr->GetLocalWeights(i));
		}

		blitz::TinyVector<double, Rank> pos;
		typename blitz::Array<cplx, Rank>::iterator it = data.begin();
		double weight = 1;
		double expValue = 0;
		for (int linearCount=0; linearCount<data.size(); linearCount++)
		{
			weight = 1;
			for (int i=0; i<Rank; i++)
			{
				pos(i) = grid(i)(it.position()(i));
				weight *= weights(i)(it.position()(i));
			}
			
			expValue += weight * std::real((*it) * conj(*it) * Potential.GetPotentialValue(pos));
			it++;
		}
		return expValue;
	}
};

#endif


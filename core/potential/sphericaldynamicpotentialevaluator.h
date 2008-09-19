#ifndef SPHERICALDYNAMICPOTENTIALEVALUATOR_H
#define SPHERICALDYNAMICPOTENTIALEVALUATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/combinedrepresentation.h"
#include "../representation/compressedrepresentation.h"
#include "staticpotential.h"
#include "potentialbase.h"
#include "potentialaction.h"

template<class DestType, class SourceType>
shared_ptr<DestType> pyprop_cast(shared_ptr<SourceType> sourcePtr)
{
#ifdef BZ_DEBUG
	shared_ptr<DestType> destPtr = dynamic_pointer_cast<DestType>(sourcePtr);
	if (destPtr == 0)
	{
		cout 
			<< "WARNING: Could not cast " << typeid(*sourcePtr).name() 
			<< " to " << typeid(DestType).name() << " safely" << endl;
		destPtr = boost::static_pointer_cast<DestType>(sourcePtr);
	}
	return destPtr;
#else
	return boost::static_pointer_cast<DestType>(sourcePtr);
#endif
}

// The class inherits from DynamicPotentialClass and ActionClass
// DynamicPotentialClass defines the potential
// ActionClass defines the how the potential is to be applied to the data
template<class ActionClass, int Rank>
class SphericalDynamicPotentialEvaluatorBase : public ActionClass
{
public:
	/** Updates propagates the wavefunction a step with this dynamic potential **/
	template<class DynamicPotentialClass>
	void IterateAction(DynamicPotentialClass &potential, const Wavefunction<Rank> &psi, blitz::Array<cplx, Rank> updateData, const cplx &timeStep, const double &curTime)
	{
		typedef CombinedRepresentation<Rank> CombRepr;		

		//Set up PotentialClass
		potential.CurTime = curTime;
		potential.TimeStep = timeStep;
		potential.CurTimeUpdated();

		//Get representations
		typename CombRepr::Ptr repr = pyprop_cast< CombRepr >(psi.GetRepresentation());
		BZPRECONDITION(repr != 0);
		typename CompressedRepresentation::Ptr angularRepr = pyprop_cast< CompressedRepresentation >(repr->GetRepresentation(Rank-1));
		BZPRECONDITION(angularRepr != 0);

		blitz::TinyVector< blitz::Array<double, 1>, Rank-1 > grid;
		for (int i=0; i<Rank-1; i++)
		{
			grid(i).reference(repr->GetLocalGrid(i));
		}
		
		blitz::Array<double, 2> omegaGrid;
		omegaGrid.reference(angularRepr->GetLocalExpandedGrid());

		//Check that omegaGrid refers to 2-d expanded grid
		BZPRECONDITION(omegaGrid.extent(1) == 2);
		
		//postition is size rank+1 since <r,(l,m)> --> <r,l,m> = pos
		blitz::TinyVector<double, Rank+1> pos;
		typename blitz::Array<cplx, Rank>::iterator it = updateData.begin();
		for (int linearCount=0; linearCount<updateData.size(); linearCount++)
		{
			for (int i=0; i<Rank-1; i++)
			{
				pos(i) = grid(i)(it.position()(i));
			}
			
			// second last rank - theta or l
			int omegaIndex = it.position()(Rank-1);
			pos(Rank-1) = omegaGrid(omegaIndex, 0);
							
			// last rank - phi or m
			pos(Rank) = omegaGrid(omegaIndex, 1);

			// Uses the function from the inherited classes
			ApplyAction(it, potential.GetPotentialValue(pos), timeStep );
			it++;
		}
	}
};

template <class PotentialClass, int Rank>
class SphericalDynamicPotentialEvaluator
{
private:
	typedef ApplyPotentialClass<Rank> ApplyClass;
	typedef MultiplyPotentialClass<Rank> MultiplyClass;
	typedef UpdatePotentialClass<Rank> UpdateClass;
	typedef GetPotentialClass<Rank> GetClass;
	
	PotentialClass Potential;
	SphericalDynamicPotentialEvaluatorBase< ApplyClass, Rank> Apply;       //Calculates exp(V)psi -> psi
	SphericalDynamicPotentialEvaluatorBase< MultiplyClass, Rank> Multiply; //Calculates V psi -> psi
	SphericalDynamicPotentialEvaluatorBase< UpdateClass, Rank> Update;     //Updates a static potential	with exp(V)
	SphericalDynamicPotentialEvaluatorBase< GetClass, Rank> Get;	       //Returns V

public:
	void ApplyConfigSection(const ConfigSection &config)
	{
		Potential.ApplyConfigSection(config);
	}

	/** Propagates the wavefunction a step with this dynamic potential **/
	void ApplyPotential(Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		Apply.IterateAction(Potential, psi, psi.GetData(), timeStep, curTime);
	}

	/** Multiplies the wavefunction with this potential. Useful for eigenproblems and expectation values **/
	void MultiplyPotential(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi, const cplx &timeStep, const double &curTime)
	{
		Multiply.DestIterator = dstPsi.GetData().begin();
		Multiply.IterateAction(Potential, srcPsi, srcPsi.GetData(), timeStep, curTime);
	}

	void UpdateStaticPotential(StaticPotential<Rank> &potential, const Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime, typename StaticPotential<Rank>::StorageModel storage)
	{
		blitz::Array<cplx, Rank> potentialData( potential.GetPotentialData() );
		if (storage == StaticPotential<Rank>::StorageExpValue)
		{
			Update.IterateAction(Potential, psi, potentialData, timeStep, curTime);
		}
		else
		{
			Get.IterateAction(Potential, psi, potentialData, timeStep, curTime);
		}
	}

	/** 
   	  * Returns an array containing this dynamic potential at the given time. 
	  * WARNING: This method allocates and returns a new blitz array in the same size as the wavefunction
	  *
	  * In a multiproc environment this returns the local portion of the potential corresponding to 
	  * the current distribution of the wavefunction.
	  */
	blitz::Array<cplx, Rank> GetPotential(const Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		blitz::Array<cplx, Rank> potentialData(psi.Data.shape());
		Get.IterateAction(Potential, psi, potentialData, timeStep, curTime);
		return potentialData;
	}

	double CalculateExpectationValue(const Wavefunction<Rank> &psi, const cplx &timeStep, const double curTime)
	{
		typedef CombinedRepresentation<Rank> CombRepr;		
	
		blitz::Array<cplx, Rank> data(psi.Data);
		
		//Set up PotentialClass
		Potential.CurTime = curTime;
		Potential.TimeStep = timeStep;
		Potential.CurTimeUpdated();

		//Get representations
		typename CombRepr::Ptr repr = dynamic_pointer_cast< CombRepr >(psi.GetRepresentation());
		typename CompressedRepresentation::Ptr angularRepr = dynamic_pointer_cast< CompressedRepresentation >(repr->GetRepresentation(Rank-1));
		BZPRECONDITION(angularRepr != 0);

		blitz::TinyVector< blitz::Array<double, 1>, Rank > grid;
		for (int i=0; i<Rank; i++)
		{
			grid(i).reference(repr->GetLocalGrid(i));
		}
		blitz::Array<double, 2> omegaGrid(angularRepr->GetLocalExpandedGrid());
		BZPRECONDITION(omegaGrid.extent(1) == 2);

		blitz::TinyVector< blitz::Array<double, 1>, Rank > weights;
		for (int i=0; i<Rank; i++)
		{
			weights(i).reference(repr->GetLocalWeights(i));
		}

		blitz::TinyVector<double, Rank+1> pos;
		typename blitz::Array<cplx, Rank>::iterator it = data.begin();
		double weight = 1;
		double expValue = 0;
		for (int linearCount=0; linearCount<data.size(); linearCount++)
		{
			weight = 1;
			for (int i=0; i<Rank-1; i++)
			{
				pos(i) = grid(i)(it.position()(i));
				weight *= weights(i)(it.position()(i));
			}
			weight *= weights(Rank-1)(it.position()(Rank-1));

			// second last rank - theta
			int omegaIndex = it.position()(Rank-1);
			pos(Rank-1) = omegaGrid(omegaIndex, 0);
							
			// last rank - phi
			pos(Rank) = omegaGrid(omegaIndex, 1);		
			
			expValue += weight * real((*it) * conj(*it) * Potential.GetPotentialValue(pos));
			it++;
		}
		return expValue;
	}

};

#endif


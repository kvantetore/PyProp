#ifndef SPHERICALDYNAMICPOTENTIALEVALUATOR_H
#define SPHERICALDYNAMICPOTENTIALEVALUATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/sphericalrepresentation3d.h"
#include "staticpotential.h"

// Action classes define the operations to be performed on the data

// Apply de Potential directly to the data
template <int Rank>
class ApplyPotentialClass
{
public:
	typedef typename blitz::Array<cplx,Rank>::iterator iterator;
		
	void ApplyAction(iterator &it, const cplx &potentialValue, const cplx &timeStep)
	{
		*it *= exp( - I * potentialValue * timeStep );
	}
};

// Update the Potential data with the potential to be applied
template <int Rank>
class UpdatePotentialClass
{
public:
	typedef typename blitz::Array<cplx,Rank>::iterator iterator;
	
	void ApplyAction(iterator &it, const cplx &potentialValue, const cplx &timeStep)
	{
		*it = exp( - I * potentialValue * timeStep );
	}
};

// Get the potential values in the data
template <int Rank>
class GetPotentialClass
{
public:
	typedef typename blitz::Array<cplx,Rank>::iterator iterator;
	
	void ApplyAction(iterator &it, const cplx &potentialValue, const cplx &timeStep)
	{
		*it = potentialValue;
	}

};



// The class inherits from DynamicPotentialClass and ActionClass
// DynamicPotentialClass defines the potential
// ActionClass defines the how the potential is to be applied to the data
template<class DynamicPotentialClass, class ActionClass, int Rank>
class SphericalDynamicPotentialEvaluatorBase : public DynamicPotentialClass, public ActionClass
{
public:
	/** Updates propagates the wavefunction a step with this dynamic potential **/
	void IterateAction(const Wavefunction<Rank> &psi, blitz::Array<cplx, Rank> updateData, const cplx &timeStep, const double &curTime)
	{
		//Set up PotentialClass
		this->CurTime = curTime;
		this->TimeStep = timeStep;

		//Get representations
		SphericalRepresentation3D *repr = static_cast<SphericalRepresentation3D*>(&psi.GetRepresentation());

		blitz::Array<double, 1> radialGrid;
		radialGrid.reference(repr->GetLocalGrid(psi, 0));

		blitz::Array<double, 2> omegaGrid;
		omegaGrid.reference(repr->GetLocalAngularGrid(psi));
		
		//postition is size rank+1 since <r,(l,m)> --> <r,l,m> = pos
		blitz::TinyVector<double, Rank+1> pos;
		typename blitz::Array<cplx, Rank>::iterator it = updateData.begin();
		for (int linearCount=0; linearCount<updateData.size(); linearCount++)
		{
			// rank 0 - radius
			pos(0) = radialGrid(it.position()(0));
			
			// rank 1 - theta
			int omegaIndex = it.position()(1);
			pos(1) = omegaGrid(omegaIndex, 0);
							
			// rank 2 - phi
			pos(2) = omegaGrid(omegaIndex, 1);

			// Uses the function from the inherited classes
			ApplyAction(it, GetPotentialValue(pos), timeStep );
			it++;
		}
	}
};

template <class PotentialClass, int Rank>
class SphericalDynamicPotentialEvaluator
{
public:
	typedef ApplyPotentialClass<Rank> ApplyClass;
	typedef UpdatePotentialClass<Rank> UpdateClass;
	typedef GetPotentialClass<Rank> GetClass;
		
	SphericalDynamicPotentialEvaluatorBase< PotentialClass, ApplyClass, Rank> Apply;
	SphericalDynamicPotentialEvaluatorBase< PotentialClass, UpdateClass, Rank> Update;
	SphericalDynamicPotentialEvaluatorBase< PotentialClass, GetClass, Rank> Get;

	
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


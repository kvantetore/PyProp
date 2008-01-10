#ifndef RANKONEPOTENTIALEVALUATOR_H
#define RANKONEPOTENTIALEVALUATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/representation.h"
#include "../utility/blitztricks.h"
#include "potentialbase.h"

/*
 * A one dimensional potential evaluator for a Rank-dimensional wavefunction.
 *
 * The potential is only dependent on one variable: along rank "potential_rank"
 * read from configSection
 *
 * The varying coordinate will be given to PotentialClass::GetPotentialValue(pos) 
 * as pos(0).
 *
 * PotentialClass may implement bool IsTimeDependent(), which returns whether the
 * potential needs to be updated every timestep. 
 */
template<class PotentialClass, int Rank>
class RankOnePotentialEvaluator
{
private:
	PotentialClass Potential;
	int PotentialRank;
	
	cplx CurTimeStep;
	blitz::Array<cplx, 1> PotentialBuffer;
	blitz::Array<cplx, 1> ExpPotentialBuffer;

	/*
	 * Setup the potential for this timestep. If required this function will update PotentialBuffer
	 * to correspond to V(x, t), if it is not timedependent this will only be done the first time the
	 * function is called
	 */
	void SetupPotentialBuffer(const Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		//Set up PotentialClass
		Potential.CurTime = curTime;
		Potential.TimeStep = timeStep;
		Potential.CurTimeUpdated();

		bool isFirstTime = (PotentialBuffer.size() == 0);
		bool isTimeDependent = Potential.IsTimeDependent();
	
		//Get grid
		typename Representation<Rank>::Ptr repr = psi.GetRepresentation();
		blitz::Array<double, 1> grid;
		grid.reference(repr->GetLocalGrid(PotentialRank));

		//Allocate if this is the first time
		if (isFirstTime)
		{
			PotentialBuffer.resize( grid.size() );
		}

		//Update Potential if it might have changed
		if (isFirstTime || isTimeDependent)
		{
			blitz::TinyVector<double, 1> pos;
			for (int i=0; i<PotentialBuffer.size(); i++)
			{
				pos(0) = grid(i);
				PotentialBuffer(i) = Potential.GetPotentialValue(pos);
			}
		}
	}

	/* 
	 * Set up the exponentiated potential exp(- i dt V)
	 *
	 * This way we eliminate a lot of exp()-calls
	 */
	void SetupExpPotentialBuffer(const Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		bool isFirstTime = (ExpPotentialBuffer.size() == 0);
		bool isTimeDependent = Potential.IsTimeDependent();
		bool timestepChanged = (timeStep != CurTimeStep);

		//Make sure we have the potential
		SetupPotentialBuffer(psi, timeStep, curTime);
		
		//Allocate if this is the first time
		if (isFirstTime)
		{
			ExpPotentialBuffer.resize( PotentialBuffer.size() );
		}

		//Update ExpPotential if it might be changed
		if (isFirstTime || isTimeDependent || timestepChanged)
		{
			CurTimeStep = timeStep;

			for (int i=0; i<ExpPotentialBuffer.size(); i++)
			{
				ExpPotentialBuffer(i) = exp(- I * PotentialBuffer(i) * timeStep);
			}
		}
	}


public:
	
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("potential_rank", PotentialRank);
		Potential.ApplyConfigSection(config);
	}

	/** Propagates the wavefunction a step with this dynamic potential **/
	void ApplyPotential(Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		//Set up potentialBuffer
		SetupExpPotentialBuffer(psi, timeStep, curTime);

		//Map the wavefunction to a rank 3 array:
		//  The first rank is all ranks greater (with greater stride) than PotentialRank
		//  The second rank is PotentialRank
		//  The third rank is all ranks less (with less stride) than PotentialRank
		blitz::Array<cplx, Rank> updateData = psi.GetData();
		blitz::Array<cplx, 3> updateData3D = MapToRank3(updateData, PotentialRank, 1);

		//We assume the data is contiguous and stored in row-major ordering (C-style)
		typename blitz::Array<cplx, 3>::iterator it = updateData3D.begin();
		for (int i=0; i<updateData3D.extent(0); i++)
		{
			for (int j=0; j<updateData3D.extent(1); j++)
			{
				cplx pot = ExpPotentialBuffer(j);
				for(int k=0; k<updateData3D.extent(2); k++)
				{
					*it *= pot;
					it++; 
				}
			}
		}
	}

	/** Multiplies the wavefunction with this potential. Useful for eigenproblems and expectation values **/
	void MultiplyPotential(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi, const cplx &timeStep, const double &curTime)
	{
		//Set up potentialBuffer
		SetupPotentialBuffer(srcPsi, timeStep, curTime);

		//Map the wavefunction to a rank 3 array:
		blitz::Array<cplx, Rank> srcData = srcPsi.GetData();
		blitz::Array<cplx, Rank> dstData = dstPsi.GetData();
		blitz::Array<cplx, 3> srcData3D = MapToRank3(srcData, PotentialRank, 1);
		blitz::Array<cplx, 3> dstData3D = MapToRank3(dstData, PotentialRank, 1);

		//We assume the data is contiguous and stored in row-major ordering (C-style)
		typename blitz::Array<cplx, 3>::iterator srcIt = srcData3D.begin();
		typename blitz::Array<cplx, 3>::iterator dstIt = dstData3D.begin();

		for (int i=0; i<srcData3D.extent(0); i++)
		{
			for (int j=0; j<srcData3D.extent(1); j++)
			{
				cplx pot = PotentialBuffer(j);
				for(int k=0; k<srcData3D.extent(2); k++)
				{
					*dstIt += pot * (*srcIt);
					dstIt++;;
					srcIt++; 
				}
			}
		}
	}

	/** 
   	  * Returns an array containing this dynamic potential at the given time. 
	  *
	  * In a multiproc environment this returns the local portion of the potential corresponding to 
	  * the current distribution of the wavefunction.
	  *
	  * Returns a 1D array corresponding to the potential values in the rank PotentialRank
	  */
	blitz::Array<cplx, 1> GetPotential(const Wavefunction<Rank> &psi, const cplx &timeStep, const double &curTime)
	{
		SetupPotentialBuffer(psi, timeStep, curTime);
		return PotentialBuffer;
	}

	cplx CalculateExpectationValue(const Wavefunction<Rank> &psi, const cplx &timeStep, const double curTime)
	{
		blitz::Array<cplx, Rank> data(psi.Data);
			
		//Set up PotentialBuffer
		SetupPotentialBuffer(psi, timeStep, curTime);

		//Get weights
		typename Representation<Rank>::Ptr repr = psi.GetRepresentation();
		blitz::TinyVector< blitz::Array<double, 1>, Rank > weights;
		for (int i=0; i<Rank; i++)
		{
			weights(i).reference(repr->GetLocalWeights(i));
		}

		blitz::TinyVector<double, Rank> pos;
		typename blitz::Array<cplx, Rank>::iterator it = data.begin();
		double weight = 1;
		cplx expValue = 0;
		for (int linearCount=0; linearCount<data.size(); linearCount++)
		{
			weight = 1;
			for (int i=0; i<Rank; i++)
			{
				weight *= weights(i)(it.position()(i));
			}
			
			expValue += weight * (*it) * conj(*it) * PotentialBuffer(it.position()(PotentialRank));
			it++;
		}
		return expValue;
	}
};


#endif


#ifndef DYNAMICPOTENTIALEVALUATOR_H
#define DYNAMICPOTENTIALEVALUATOR_H

template<class DynamicPotentialClass, int Rank>
class DynamicPotentialEvaluator : public DynamicPotentialClass
{
public:
	DynamicPotentialEvaluator() {}
	virtual ~DynamicPotentialEvaluator() {}

	/** Updates propagates the wavefunction a step with this dynamic potential **/
	dbl ApplyPotential(blitz::Array<cplx, Rank> &psi, cplx timeStep, dbl curTime)
	{
		//Set up PotentialClass
		this->CurTime = curTime;
		this->TimeStep = timeStep;
	
		//Set up grid
		blitz::TinyVector< blitz::Array<dbl, 1>, Rank> grid;
		for (int curRank = 0; curRank<Rank; curRank++)
		{
			grid(curRank).resize(psi.extent(curRank));
			grid(curRank) = blitz::tensor::i;
		}
	
		blitz::TinyVector<dbl, Rank> pos;
		blitz::TinyVector<int, Rank> indexPos;
		typename blitz::Array<cplx, Rank>::iterator it = psi.begin();
		
		HL::Timer timer;
		timer.start();


		for (int linearCount=0; linearCount<psi.size(); linearCount++)
		{
			for (int curRank=0; curRank<Rank; curRank++)
			{
				pos(curRank) = grid(curRank)(it.position()(curRank));
			}
			cplx potValue = GetPotentialValue(pos);
			*it *= exp( - I * potValue * timeStep );

			it++;
		}


/*		//Standard C++
		int linearCount = 0;
		cplx* data = psi.data();
		dbl potentialValue = 0.0;
		for (int i=0;i<psi.extent(0);i++)
		{
			for (int j=0;j<psi.extent(1);j++)
			{
				pos(0) = grid(0)(j);
				pos(1) = grid(1)(j);
	
				potentialValue = GetPotentialValue(pos);
				data[linearCount] *= exp(-I * timeStep * potentialValue) ;
				linearCount++;
			}
		}

*/		
/*	
		//Blitz tensor notation
		blitz::firstIndex i;
		blitz::secondIndex j;
		psi = psi(i,j) * exp( -I * timeStep * (sqr(grid(0)(i)) + sqr(grid(1)(j))));
*/

		
		timer.stop();
		dbl timeUsed = (dbl)timer;
		dbl gps = (dbl)psi.size() / timeUsed;
		return gps;
	}
};

/* 
Dynamic potential for evaluation of the kinetic energy potential for CartesianFFTEvaluator
*/
template<int Rank>
class CartesianKineticEnergyPotential
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	dbl CurTime;

	dbl var;
	
	void ApplyConfigSection(const ConfigSection &config)
	{
		var = 100.3;
	}

	dbl GetPotentialValue(const blitz::TinyVector<dbl, Rank> &momentum)
	{
		dbl kineticPotential = 0.0;

		for (int i=0;i<Rank;i++)
		{
			kineticPotential += (momentum(i) * momentum(i)) ;
		}
		return kineticPotential/2.0;
	}
};

#endif


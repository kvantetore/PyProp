
#include <boost/python.hpp>
#include "../common.h"
#include "../wavefunction.h"
#include "../representation/representation.h"
#include "../potential/staticpotential.h"

using namespace boost::python;

template<int Rank>
void SetWavefunctionFromGridFunction(Wavefunction<Rank> &psi, object& function, object &conf)
{
	//Set up grid
	blitz::TinyVector< blitz::Array<double, 1>, Rank> grid;
	typename Representation<Rank>::Ptr repr = psi.GetRepresentation();
	for (int curRank = 0; curRank<Rank; curRank++)
	{
		grid(curRank).reference(repr->GetLocalGrid(curRank));
	}
	
	//blitz::TinyVector<double, Rank> pos;
	list pos;
	for (int curRank=0; curRank<Rank; curRank++)
	{
		pos.append(0.0);
	}
	blitz::TinyVector<int, Rank> indexPos;
	typename blitz::Array<cplx, Rank>::iterator it = psi.Data.begin();
	for (int linearCount=0; linearCount<psi.Data.size(); linearCount++)
	{
		for (int curRank=0; curRank<Rank; curRank++)
		{
			pos[curRank] = grid(curRank)(it.position()(curRank));
		}
		
		//cplx value = extract<cplx>(function(conf, pos));
		cplx value = extract<cplx>(function(conf,pos));
		*it = value;
		
		it++;
	}

}

template<int Rank>
void SetPotentialFromGridFunction(StaticPotential<Rank> &potential, cplx timeStep, Wavefunction<Rank> &psi, Representation<Rank> &repr, object function, object conf)
{
	//Set up grid
	blitz::TinyVector< blitz::Array<double, 1>, Rank> grid;
	for (int curRank = 0; curRank<Rank; curRank++)
	{
		grid(curRank).reference(repr.GetLocalGrid(curRank));
	}
	
	list pos;
	for (int curRank=0; curRank<Rank; curRank++)
	{
		pos.append(0.0);
	}
	blitz::TinyVector<int, Rank> indexPos;
	typename blitz::Array<cplx, Rank>::iterator it = potential.GetPotentialData().begin();
	for (int linearCount=0; linearCount<potential.GetPotentialData().size(); linearCount++)
	{
		for (int curRank=0; curRank<Rank; curRank++)
		{
			pos[curRank] = grid(curRank)(it.position()(curRank));
		}
		
		cplx value = extract<cplx>(function(conf, pos));
		*it = exp( - I * value * timeStep );
		
		it++;
	}
}


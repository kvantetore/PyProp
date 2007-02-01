#ifndef CARTESIANKINETICPOTENTIAL_H
#define CARTESIANKINETICPOTENTIAL_H

#include "../common.h"
#include "../wavefunction.h"
#include "../representation/cartesianrepresentation.h"

template<int Rank>
class CartesianKineticPotential
{
public:
	blitz::Array<cplx, Rank> ExpPotential;

	CartesianKineticPotential(CartesianRepresentation<Rank> &repr, double dt)
	{
		//Allocate Potential Array
		blitz::TinyVector<int, Rank> shape = repr.GetInitialShape();
		ExpPotential.resize(shape);
		
		//Gather grid arrays
		blitz::TinyVector< blitz::Array<double, 1>, Rank> grid;
		for (int i=0; i<Rank; i++)
		{
			grid(i).reference(repr.Range(i).GetGrid());
		}
		
		CreatePotential(grid, dt);
	}
	
	void CreatePotential(blitz::TinyVector< blitz::Array<double, 1>, Rank> &grid, double dt);
	
	void ApplyPotential(Wavefunction<Rank> &psi)
	{
		psi.Data = psi.Data * ExpPotential;
	}
};

template<int Rank>
void CartesianKineticPotential<Rank>::CreatePotential(blitz::TinyVector< blitz::Array<double, 1>, Rank> &grid, double dt)
{
	std::cerr << "CartesianKinetPotential not implemented for Rank " << Rank << " > 3" << std::endl;
}

template<>
void CartesianKineticPotential<1>::CreatePotential(blitz::TinyVector< blitz::Array<double, 1>, 1> &grid, double dt)
{
	blitz::Array<double, 1> k(grid(0));
	ExpPotential = exp(- I * k * k * dt);
}

template<>
void CartesianKineticPotential<2>::CreatePotential(blitz::TinyVector< blitz::Array<double, 1>, 2> &grid, double dt)
{
	using namespace blitz::tensor;
	
	blitz::Array<double, 1> k0(grid(0));
	blitz::Array<double, 1> k1(grid(1));
	
	ExpPotential = exp(- I * ((k0(i)* k0(i)) + (k1(j) * k1(j))) * dt);
}

template<>
void CartesianKineticPotential<3>::CreatePotential(blitz::TinyVector< blitz::Array<double, 1>, 3> &grid, double dt)
{
	using namespace blitz::tensor;
	
	blitz::Array<double, 1> k0(grid(0));
	blitz::Array<double, 1> k1(grid(1));
	blitz::Array<double, 1> k2(grid(2));
	
	ExpPotential = exp(- I * ((k0(i)* k0(i)) + (k1(j) * k1(j)) + (k2(k) * k2(k))) * dt);
}


#endif


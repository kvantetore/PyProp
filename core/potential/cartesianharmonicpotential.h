#ifndef CARTESIANHARMONICPOTENTIAL
#define CARTESIANHARMONICPOTENTIAL

#include "../representation/cartesianrepresentation.h"
#include "../wavefunction.h"

template<int Rank>
class CartesianHarmonicPotential
{
private:
	double Strength;

	void inline ApplyDynamicPotential(Wavefunction<Rank> &wave, blitz::TinyVector< blitz::Array<double, 1>, Rank> &grid, double t, double dt);
	
public:

	CartesianHarmonicPotential(double strength) :
		Strength(strength)
	{ }

	void ApplyDynamicPotential(Wavefunction<Rank> &wave, double t, double dt)
	{
		CartesianRepresentation<Rank> *repr = dynamic_cast< CartesianRepresentation<Rank>* >(&wave.GetRepresentation());

		//Gather grid arrays
		blitz::TinyVector< blitz::Array<double, 1>, Rank> grid;
		for (int i=0; i<Rank; i++)
		{
			grid(i).reference(repr->GetLocalGrid(wave, i));
		}
		
		ApplyDynamicPotential(wave, grid, t, dt);
	}

};

template<int Rank> 
void inline CartesianHarmonicPotential<Rank>::ApplyDynamicPotential(Wavefunction<Rank> &wave, blitz::TinyVector< blitz::Array<double, 1>, Rank> &grid, double t, double dt)
{
	std::cerr << "Cartesian Harmonic Potential not implemented for rank " << Rank << std::endl;
}

template<> 
void inline CartesianHarmonicPotential<1>::ApplyDynamicPotential(Wavefunction<1> &wave, blitz::TinyVector< blitz::Array<double, 1>, 1> &grid, double t, double dt)
{
	using namespace blitz;
	using namespace blitz::tensor;
	Array<double, 1> x(grid(firstDim));
	
	wave.Data = wave.Data * exp(- I * dt * Strength * sqr(x) );
}

template<> 
void inline CartesianHarmonicPotential<2>::ApplyDynamicPotential(Wavefunction<2> &wave, blitz::TinyVector< blitz::Array<double, 1>, 2> &grid, double t, double dt)
{
	using namespace blitz;
	using namespace blitz::tensor;
	Array<double, 1> x(grid(firstDim));
	Array<double, 1> y(grid(firstDim));
	
	wave.Data = wave.Data * exp(- I * dt * Strength * (sqr(x(i)) + sqr(y(j))) );
}

template<> 
void inline CartesianHarmonicPotential<3>::ApplyDynamicPotential(Wavefunction<3> &wave, blitz::TinyVector< blitz::Array<double, 1>, 3> &grid, double t, double dt)
{
	using namespace blitz;
	using namespace blitz::tensor;
	Array<double, 1> x(grid(firstDim));
	Array<double, 1> y(grid(secondDim));
	Array<double, 1> z(grid(thirdDim));
	
	wave.Data = wave.Data * exp(- I * dt * Strength * (sqr(x(i)) + sqr(y(j)) + sqr(z(k))) );
}


#endif

#ifndef ORTHOPOLPOTENTIALS_H
#define ORTHOPOLPOTENTIALS_H

#include "../common.h"
#include "../wavefunction.h"
#include "dynamicpotentialevaluator.h"

template<int Rank>
class AddedHermitePotential
{
public:
	double CurTime;
	cplx TimeStep;
	double Alpha;
	int RadialRank;
	int hypersphericalDimensions;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("alpha", Alpha);
		config.Get("radial_rank", RadialRank);
		config.Get("hyperspherical_dimensions", hypersphericalDimensions);
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{	
		double r = pos(RadialRank);
		double x = Alpha * r;
		//int N = hypersphericalDimensions;
		//return sqr(Alpha) * ( 1.0 + (N-1)*(N-3)/(sqr(r) * 4.0) - sqr(x) ) / 2; 
		return sqr(Alpha) * ( 1.0 - sqr(x) ) / 2.0; 
	}
};

template<int Rank>
class AddedLaguerrePotential
{
public:
	double CurTime;
	cplx TimeStep;
	double Alpha;
	int RadialRank;
	int hypersphericalDimensions;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("alpha", Alpha);
		config.Get("radial_rank", RadialRank);
		config.Get("hyperspherical_dimensions", hypersphericalDimensions);
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{	
		double r = pos(RadialRank);
		double x = sqr( Alpha * r ) / 4.0;
		//int N = hypersphericalDimensions;
		//return sqr(Alpha) * ( N - x ) / 8.0; 
		return sqr(Alpha) * ( x + 1.0 ) / 8.0; 
	}
};

#endif


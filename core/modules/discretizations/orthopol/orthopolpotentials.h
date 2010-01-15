#ifndef ORTHOPOLPOTENTIALS_H
#define ORTHOPOLPOTENTIALS_H

#include "../common.h"
#include "../wavefunction.h"
#include "rankonepotentialevaluator.h"

template<int Rank>
class AddedHermitePotential : public PotentialBase<Rank>
{
public:
	double CurTime;
	cplx TimeStep;
	double Scaling;
	int hypersphericalDimensions;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("scaling", Scaling);
		config.Get("hyperspherical_dimensions", hypersphericalDimensions);
	}

	bool IsTimeDependent()
	{
		return false;
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{	
		double r = pos(0);
		double x = Scaling * r;
		//int N = hypersphericalDimensions;
		//return sqr(Scaling) * ( 1.0 + (N-1)*(N-3)/(sqr(r) * 4.0) - sqr(x) ) / 2; 
		return sqr(Scaling) * ( 1.0 - sqr(x) ) / 2.0; 
	}
};

template<int Rank>
class AddedLaguerrePotential : public PotentialBase<Rank>
{
public:
	double CurTime;
	cplx TimeStep;
	double Scaling;
	int hypersphericalDimensions;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("scaling", Scaling);
		config.Get("hyperspherical_dimensions", hypersphericalDimensions);
	}
	
	bool IsTimeDependent()
	{
		return false;
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{	
		/*
		double r = pos(0);
		double x = sqr( Scaling * r ) / 4.0;
		//int N = hypersphericalDimensions;
		//return sqr(Scaling) * ( N - x ) / 8.0; 
		return sqr(Scaling) * ( x + 1.0 ) / 8.0; 
		*/

		const int D = hypersphericalDimensions;
		const double x = sqr(Scaling * pos(0)) / 4.0;

		return sqr(Scaling) * (D - x) / 8.0;
	}
};

#endif


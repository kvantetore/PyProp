#ifndef DIAGONALPOTENTIALS_H
#define DIAGONALPOTENTIALS_H

#include "../common.h"
#include "../wavefunction.h"
#include "sphericaldynamicpotentialevaluator.h"

// Potential classes to be used by the SphericalDynamicPotentialEvaluator

template<int Rank>
class DiagonalAngularPotential
{

public:
	double CurTime;
	cplx TimeStep;

	void ApplyConfigSection(const ConfigSection &config)
	{
		
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{
		double r = pos(0);
		if (fabs(r) < 1e-5) {
			return 0;
		}
		double l = pos(1);
		//int m = (int)pos(2);
		return l*(l+1.0)/(2.0*r*r);
	}
};

template<int Rank>
class DiagonalRadialPotential
{
public:
	double CurTime;
	cplx TimeStep;
	
	void ApplyConfigSection(const ConfigSection &config)
	{
		
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{
		return sqr(pos(0))/2;
	
	}
};

template <int Rank>
class InitialPotential
{
public:
	double CurTime;
	cplx TimeStep;

	void ApplyConfigSection(const ConfigSection &config)
	{
		
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(0);
		if (fabs(r) < 1e-10) {
			return 0;
		}
		return -1/fabs(r);
	}
};


#endif


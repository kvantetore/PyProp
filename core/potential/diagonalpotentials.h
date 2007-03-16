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
	double Mass;
	int RadialRank;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", Mass);
		config.Get("radial_rank", RadialRank);
		cout << "DiagonalAngularPotential: mass = " << Mass << ", radial_rank = " << RadialRank << endl;
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{
		double r = pos(RadialRank);
		if (fabs(r) < 1e-5) {
			return 0;
		}
		double l = pos(Rank-2);
		//int m = (int)pos(Rank-1);
		return l*(l+1.0)/(2.0*Mass*r*r);
	}
};

template<int Rank>
class DiagonalRadialPotential
{
public:
	double CurTime;
	cplx TimeStep;
	double Mass;
	int RadialRank;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", Mass);
		config.Get("radial_rank", RadialRank);
		cout << "DiagonalRadialPotential: mass = " << Mass << ", radial_rank = " << RadialRank << endl;
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos )
	{
		return sqr(pos(RadialRank))/(2*Mass);
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
		throw std::runtime_error("Invalid potential");
	}

	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = pos(Rank-3);
		if (fabs(r) < 1e-10) {
			return 0;
		}
		return -1/fabs(r);
	}
};


#endif


#ifndef RADIALTRANSFORM_H
#define RADIALTRANSFORM_H

#include <core/common.h>
#include <core/wavefunction.h>

template <int Rank>
class RadialTransform{
public:
	RadialTransform() {}
	~RadialTransform() {}

	//Standard transform interface
	void TransformRank(Wavefunction<Rank> &psi, int rank, int direction);
	void ForwardTransform(Wavefunction<Rank> &psi, int rank);
	void InverseTransform(Wavefunction<Rank> &psi, int rank);
};

#endif

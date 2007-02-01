#ifndef CARTESIANFOURIERTRANSFORM_H
#define CARTESIANFOURIERTRANSFORM_H

#include "../common.h"
#include "../wavefunction.h"
#include "fouriertransform.h"
#include "../representation/cartesianrepresentation.h"

template<int Rank> class CartesianRepresentation;

/**
*/
template<int Rank>
class CartesianFourierTransform
{
public:
	CartesianFourierTransform() 		
	{
	}

	/*
	Simple interface 
	These functions normalize the wavefunction after the inverse transform
	*/
	void ForwardTransform(Wavefunction<Rank> &psi);	
	void InverseTransform(Wavefunction<Rank> &psi);

	/* Advanced interface */
	void TransformExceptDistributedRank(Wavefunction<Rank> &psi, int direction);
	void TransformRank(Wavefunction<Rank> &psi, int rank, int direction);
	
	/**
	Renormalizes the wavefunction after a forward/backward transform.
	Make sure to call this function after a full forward/backward transform,
	otherwise, the wavefunction will be scaled by Mul(strides)
	*/
	void Renormalize(Wavefunction<Rank> &psi);

	/**
	 * Methods for manipulating representations
	*/
	// Create representations for a full transform 
	CartesianRepresentation<Rank> CreateFourierRepresentation(const CartesianRepresentation<Rank> &gridRepresentation);
	//Create representations for partial transforms 
	CartesianRepresentation<Rank> CreateFourierRepresentation(const CartesianRepresentation<Rank> &gridRepresentation, int rank);
}
;

#endif

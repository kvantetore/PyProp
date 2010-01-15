#ifndef CARTESIANFOURIERTRANSFORM_H
#define CARTESIANFOURIERTRANSFORM_H

#include "../common.h"
#include "../wavefunction.h"
#include "fouriertransform.h"
#include "../representation/cartesianrepresentation.h"

template<int Rank> class CartesianRepresentation;

/*
 * This class really has to different interfaces: 
 * 1) The high level interface. Calling ForwardTransform() and 
 *    InverseTransform() will fourier transform the entire wavefunction,
 *    taking care of distribution among processors etc.
 *
 * 2) The low level interace. TransformRank, etc. transformes the 
 *    wavefunction along a single rank. It will return an error if the transformed
 *    rank is currently distributed, and will not perform any redistributions.
 */
template<int Rank>
class CartesianFourierTransform
{
private:
	void FourierTransform(Wavefunction<Rank> &psi, int direction);
	
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
	typename CartesianRepresentation<Rank>::Ptr CreateFourierRepresentation(const CartesianRepresentation<Rank> &gridRepresentation);
	//Create representations for partial transforms 
	typename CartesianRepresentation<Rank>::Ptr CreateFourierRepresentation(const CartesianRepresentation<Rank> &gridRepresentation, int rank);
}
;

#endif

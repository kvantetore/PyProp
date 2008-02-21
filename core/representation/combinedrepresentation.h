#ifndef COMBINEDREPRESENTATION_H
#define COMBINEDREPRESENTATION_H

#include "../wavefunction.h"
#include "representation.h"

#include <boost/shared_ptr.hpp>

typedef Representation<1>::Ptr Representation1DPtr;

template<int Rank>
class CombinedRepresentation : public Representation<Rank>
{
private:
	//Variables
	Representation1DPtr Representations[Rank];
	blitz::Array<double, 1> LocalGrid[Rank];

public:
	typedef shared_ptr< CombinedRepresentation<Rank> > Ptr;

	//Constructors
	CombinedRepresentation();
	virtual ~CombinedRepresentation();
	CombinedRepresentation(const CombinedRepresentation<Rank> &other);
	CombinedRepresentation<Rank>& operator=(const CombinedRepresentation<Rank> &other);
	
	virtual typename Representation<Rank>::RepresentationPtr Copy();

	//Get/Set the representation of the specified rank
	Representation1DPtr GetRepresentation(int rank);
	void SetRepresentation(int rank, Representation1DPtr repr);
	
	//Implementation of the representation interface.
	virtual blitz::TinyVector<int, Rank> GetFullShape();
	virtual cplx InnerProduct(const Wavefunction<Rank> &w1, const Wavefunction<Rank> &w2);
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank);
	virtual blitz::Array<double, 1> GetLocalWeights(int rank);
	virtual void ApplyConfigSection(const ConfigSection &config);
};

#endif


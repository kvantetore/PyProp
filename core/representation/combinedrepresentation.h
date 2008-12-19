#ifndef COMBINEDREPRESENTATION_H
#define COMBINEDREPRESENTATION_H

#include "../wavefunction.h"
#include "representation.h"

#include <boost/shared_ptr.hpp>

typedef Representation<1>::Ptr Representation1DPtr;

template<int Rank>
class CombinedRepresentation : public Representation<Rank>
{
public:
	typedef blitz::Array<cplx, Rank> DataArray;

private:
	//Variables
	blitz::TinyVector<Representation1DPtr, Rank> Representations;

	cplx InnerProductImpl_Algo1(DataArray d1, DataArray d2);
	cplx InnerProductImpl_Algo2(DataArray d1, DataArray d2);

public:
	typedef shared_ptr< CombinedRepresentation<Rank> > Ptr;

	int Algorithm;

	//Constructors
	CombinedRepresentation();
	virtual ~CombinedRepresentation();
	
	virtual typename Representation<Rank>::RepresentationPtr Copy();

	//Get/Set the representation of the specified rank
	Representation1DPtr GetRepresentation(int rank);
	void SetRepresentation(int rank, Representation1DPtr repr);
	
	//Implementation of the representation interface.
	virtual blitz::TinyVector<int, Rank> GetFullShape();
	virtual cplx InnerProduct(const Wavefunction<Rank> &w1, const Wavefunction<Rank> &w2);
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank);
	virtual blitz::Array<double, 1> GetGlobalWeights(int rank);
	virtual void ApplyConfigSection(const ConfigSection &config);

	/*
	 * Overlap Matrix related stuff
	 */
	virtual OverlapMatrix::Ptr GetGlobalOverlapMatrix(int rank);
	virtual void MultiplyOverlap(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi, int rank);
	virtual void MultiplyOverlap(cplx sourceScaling, Wavefunction<Rank> &srcPsi, cplx destScaling, Wavefunction<Rank> &dstPsi, int rank);
	virtual void MultiplyOverlap(Wavefunction<Rank> &psi);
	virtual void SolveOverlap(Wavefunction<Rank> &psi);
	virtual void MultiplySqrtOverlap(bool conjugate, Wavefunction<Rank> &psi);
	virtual void SolveSqrtOverlap(bool conjugate, Wavefunction<Rank> &psi);
};

#endif


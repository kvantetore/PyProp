#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "common.h"
#include <vector>

template<int Rank> class Representation;

template<int Rank> 
class Wavefunction
{
public:
	typedef boost::shared_ptr< Representation<Rank> > RepresentationPtr;
	typedef blitz::Array<cplx, Rank> DataArray;
	typedef blitz::TinyVector<int, Rank> IndexVector;
	typedef boost::shared_ptr< DataArray > DataArrayPtr;
	typedef std::vector<DataArrayPtr> DataArrayVector;
		
private:
	RepresentationPtr Repr;
	DataArrayVector WavefunctionData;
	int ActiveBufferName;
			
public:
	blitz::Array<cplx, Rank> Data;
	
	Wavefunction() 
	{
	}
	
	void SetRepresentation(RepresentationPtr repr)
	{
		Repr = repr;
	}

	Representation<Rank>& GetRepresentation() const
	{
		return *Repr;
	}
	
	blitz::Array<cplx, Rank> GetData()
	{
		return Data;
	}

	const blitz::Array<cplx, Rank> GetData() const
	{
		return Data;
	}

	/* Returns the number of bytes currently allocated by this wavefunction */
	size_t GetMemoryFootprint() const;

	/* Allocates a new wavefunction buffer of size specifed by shape,
	 * and returns the name of the data buffer
	 */
	int AllocateData(blitz::TinyVector<int, Rank> shape);

	/* Frees a previously allocated data buffer specified by it name
	 */
	void FreeData(int bufferName);

	/* Returns the name of the active data buffer */
	int GetActiveBufferName() const;

	/* Set the currently active databuffer */
	int SetActiveBuffer(int bufferName);

	/* Get a reference to a buffer */
	DataArray& GetData(int bufferName);
	const DataArray& GetData(int bufferName) const;
	
	int GetRank() const
	{
		return Rank;
	}

	void AllocateData();
	
	/* Utility functions */
	double GetNorm() const;
	double GetLocalNorm() const
	{
		return GetRepresentation().InnerProduct(*this, *this).real();
	}
	
	double Normalize()
	{
		double norm = GetNorm();
		Data /= sqrt(norm);
		return norm;
	}
	
	cplx InnerProduct(Wavefunction<Rank> &psi) const
	{
		return GetRepresentation().InnerProduct(psi, *this);
	}

	/* Make a copy of the currently active data buffer only.
	 * Use this method when you want a copy of the wavefunction for
	 * Computing values throughout the propagation (i.e. autocorrelation)
	 * Use CopyDeep() if you need to propagate on the copied wavefunction.
	 */
	Wavefunction<Rank>* Copy() const;
	Wavefunction<Rank>* CopyDeep() const;
	
};


#endif


#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "common.h"
#include "utility/databuffer.h"
#include <vector>

template<int Rank> class Representation;

template<int Rank> 
class Wavefunction
{
public:
	typedef boost::shared_ptr< Wavefunction<Rank> > Ptr;
	typedef boost::shared_ptr< Representation<Rank> > RepresentationPtr;
	typedef blitz::Array<cplx, Rank> DataArray;
	typedef blitz::TinyVector<int, Rank> IndexVector;
	typedef std::vector< DataBuffer<Rank> > DataBufferList;
		
	blitz::Array<cplx, Rank> Data;
	
	Wavefunction()
	{
		ActiveBufferName = -1;
	}
	
	void SetRepresentation(RepresentationPtr repr)
	{
		Repr = repr;
	}

	RepresentationPtr GetRepresentation() const
	{
		return Repr;
	}
	
	blitz::Array<cplx, Rank>& GetData()
	{
		return Data;
	}

	const blitz::Array<cplx, Rank>& GetData() const
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

	/* Set the current buffer name to reference another data buffer*/
	void SetData(DataArray &newData);

	/* Methods for locking and unlocking data buffers */
	void LockBuffer(int bufferName);
	void UnLockBuffer(int bufferName);
	int GetAvailableDataBufferName(const blitz::TinyVector<int, Rank> &shape) const;
	bool HasAvailableBuffer(const blitz::TinyVector<int, Rank> &shape) const;

	int GetRank() const
	{
		return Rank;
	}

	void AllocateData();
	
	/* Utility functions */
	cplx InnerProduct(const Wavefunction<Rank> &psi) const;
	cplx LocalInnerProduct(const Wavefunction<Rank> &psi) const;
	double GetNorm() const;
	double GetLocalNorm() const;
	
	double Normalize()
	{
		double norm = GetNorm();
		Data /= norm;
		return norm;
	}

	void Clear()
	{
		Data = 0;
	}
	
	/* Make a copy of the currently active data buffer only.
	 * Use this method when you want a copy of the wavefunction for
	 * Computing values throughout the propagation (i.e. autocorrelation)
	 * Use CopyDeep() if you need to propagate on the copied wavefunction.
	 */
	Ptr Copy() const;
	Ptr CopyDeep() const;
	
private:
	RepresentationPtr Repr;
	DataBufferList WavefunctionData;
	int ActiveBufferName;

};


#endif


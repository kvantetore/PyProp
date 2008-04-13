#ifndef DATABUFFER_H
#define DATABUFFER_H

#include "../common.h"

/*
 * A Wavefunction maintains a list of DataBuffers. Each DataBuffer is a complex array
 * of the same rank as the wavefunction. Buffers are allocated by algorithms whenever temporary
 * buffers are needed. The algorithm then locks the DataBuffer, indicating that other algorithms should
 * not use this DataBuffer. When the algorithm has completed its use of the DataBuffer, the lock is released, 
 * and other algorithms can then use the same DataBuffer.
 */

template<int Rank>
class DataBuffer
{
public:
	typedef boost::shared_ptr<bool> BoolPtr;
	typedef blitz::Array<cplx, Rank> ArrayType;

	DataBuffer()
	{
		DataArray.resize(0);
		isAvailable = BoolPtr(new bool);
		*isAvailable = true;
	}

	DataBuffer(blitz::TinyVector<int, Rank> &shape)
	{
		DataArray.resize(shape);
		isAvailable = BoolPtr(new bool);
		*isAvailable = true;
	}

	DataBuffer(const DataBuffer<Rank> &other)
	{
		DataArray.reference(other.DataArray);
		isAvailable = other.isAvailable;
	}

	DataBuffer<Rank>& operator=(blitz::TinyVector<int, Rank> &other) 
	{
		DataArray.reference(other.DataArray);
		isAvailable = other.isAvailable;
		return *this;
	}

	bool operator==(const blitz::TinyVector<int, Rank> &shape) const
	{
		for (int i=0; i<Rank; i++)
		{
			if (DataArray.extent(i) != shape(i))
			{
				return false;
			}
		}
		return true;
	}

	bool IsAvailable() const
	{
		return *isAvailable;
	}

	void Lock()
	{
		if (!IsAvailable())
		{
			throw std::runtime_error("Array is already in use");
		}
		*isAvailable = false;
	}

	void UnLock()
	{
		*isAvailable = true;
	}

	ArrayType& GetArray()
	{
		return DataArray;
	}

	const ArrayType& GetArray() const
	{
		return DataArray;
	}


	void ResizeArray(const blitz::TinyVector<int, Rank>& shape)
	{
		DataArray.resize(shape);
	}

	void FreeArray()
	{
		DataArray.resize(0);
	}


private:
	ArrayType DataArray;
	BoolPtr isAvailable;
};


#endif


#ifndef BLITZMPI_H
#define BLITZMPI_H

#include <mpi.h>

#include "../common.h"
#include "mpitraits.h"
#include "../utility/blitzutils.h"

#include <map>
#include <vector>

using namespace std;

template<class T, int Rank>
class BlitzMPI
{
public:
	typedef blitz::Array<T, Rank> DataArray;
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::TinyVector<MPI_Datatype, Rank> DatatypeVector;
	typedef MPITraits<T> Traits;

	typedef std::pair< DataVector, DataVector > DatatypeKey;
	typedef std::pair< DatatypeKey, MPI_Datatype > DatatypeEntry;
	typedef std::vector< DatatypeEntry > DatatypeList;

private:
	static DatatypeList Datatypes;

public:
	BlitzMPI(const DataVector &shape, const DataVector &stride)
	{
		SetupDatatype(shape, stride);
	}
	
	BlitzMPI(const DataArray &shapeData)
	{
		SetupDatatype(shapeData.shape(), shapeData.stride());
	}

	~BlitzMPI()
	{
		FreeDatatype();
	}
	
	int Send(const DataArray &data, int dest, int tag, MPI_Comm comm);
	int Recv(DataArray &data, int dest, int tag, MPI_Comm comm);
	MPI_Request ISend(const DataArray &data, int dest, int tag, MPI_Comm comm);
	MPI_Request IRecv(DataArray &data, int dest, int tag, MPI_Comm comm);

	MPI_Datatype GetDatatype()
	{
		return Datatype;
	}

private:
	MPI_Datatype Datatype;

public:
	bool IsContiguous(const DataVector &shape, const DataVector &stride);
	void SetupDatatype(const DataVector &shape, const DataVector &stride);
	void FreeDatatype();
};

template<class T, int Rank> typename BlitzMPI<T, Rank>::DatatypeList BlitzMPI<T, Rank>::Datatypes;


template<class T, int Rank>
bool BlitzMPI<T, Rank>::IsContiguous(const DataVector &shape_, const DataVector &stride_)
{
	int expectedStride = 1;
	for (int i=Rank-1; i>=0; i--)
	{
		if (expectedStride != stride_(i))
		{
			return false;
		}
		expectedStride *= shape_(i);
	}
	return true;
}

template<class T, int Rank>
void BlitzMPI<T, Rank>::SetupDatatype(const DataVector &shape, const DataVector &stride)
{
	//Check if we have already created the datatype
	
	bool foundDatatype = false;
	for (typename DatatypeList::iterator it=Datatypes.begin(); it!=Datatypes.end(); it++)
	{
		DataVector curShape = (*it).first.first;
		DataVector curStride = (*it).first.second;

		if ((curShape == shape) && (curStride == stride))
		{
			foundDatatype = true;
			Datatype =(*it).second;
		}
	}

	if (!foundDatatype)
	{
		int size = 1;
		for (int i=0; i<Rank; i++)
		{
			size *= shape(i);
		}
		
		if (IsContiguous(shape, stride))
		{
			/*
			 * Update: Actually, it was a bug in OpenMPI 1.1.0 (!), it works perfectly on other MPI 
			 * implementations (including OpenMPI >= 1.1.2)
			 * 
			 *
			 * I cannot, for the love of god, make this work for array size (2, 1), stride (2, 1) when sending to 
			 * array size (2, 1), stride (1, 1) or vice versa for complex numbers (using std::complex). 
			 * 
			 * If anyone can figure out why, I would be very grateful. I suspect there are some weird alignment issues,
			 * but I am not sure. -Tore
			 */
		
			MPI_Datatype prevType = Traits::Type() ;
			int origSize = Traits::Length();
			MPI_Type_contiguous(size * origSize, prevType, &Datatype);
			MPI_Type_commit(&Datatype);
		}
		else 
		{
			MPI_Datatype prevType = Traits::Type();
			MPI_Datatype curType;
		
			//first rank
			int origSize = Traits::Length();
			MPI_Type_hvector(shape(Rank-1), origSize, stride(Rank-1) * sizeof(T), prevType, &curType);
			MPI_Type_commit(&curType);
			prevType = curType;
		
			//the other ranks
			for (int i=Rank-2; i>=0; i--)
			{
				MPI_Type_hvector(shape(i), 1, stride(i)*sizeof(T), prevType, &curType);
				MPI_Type_commit(&curType);
		
				MPI_Type_free(&prevType);
				prevType = curType;
			}
		
			//store the final datatype
			int size = 0;
			MPI_Aint extent = 0;
			MPI_Type_size(curType, &size);
			MPI_Type_extent(curType, &extent);
			Datatype = curType;
		}


		DatatypeEntry entry;
		entry.first.first = shape;
		entry.first.second = stride;
		entry.second = Datatype;
		Datatypes.push_back(entry);
	}
}
	
template<class T, int Rank>
void BlitzMPI<T, Rank>::FreeDatatype()
{
	//Dont do that...
	//MPI_Type_free(&Datatype);
}

template<class T, int Rank>
int BlitzMPI<T, Rank>::Send(const DataArray &data, int dest, int tag, MPI_Comm comm)
{
	return MPI_Send(const_cast<T*>(data.data()), 1, Datatype, dest, tag, comm);
}

template<class T, int Rank>
int BlitzMPI<T, Rank>::Recv(DataArray &data, int src, int tag, MPI_Comm comm)
{
	MPI_Status status;
	return MPI_Recv(data.data(), 1, Datatype, src, tag, comm, &status);
}

template<class T, int Rank>
MPI_Request BlitzMPI<T, Rank>::ISend(const DataArray &data, int dest, int tag, MPI_Comm comm)
{
	MPI_Request request;
	MPI_Isend(const_cast<T*>(data.data()), 1, Datatype, dest, tag, comm, &request);
	return request;
}

template<class T, int Rank>
MPI_Request BlitzMPI<T, Rank>::IRecv(DataArray &data, int src, int tag, MPI_Comm comm)
{
	MPI_Request request;
	MPI_Irecv(data.data(), 1, Datatype, src, tag, comm, &request);
	return request;
}



#endif


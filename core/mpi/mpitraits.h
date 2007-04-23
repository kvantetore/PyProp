#ifndef MPITRAITS_H
#define MPITRAITS_H

#include "../common.h"
#include <mpi.h>

/*
template<class T> 
class MPITraits
{
public:
	static MPI_Datatype Type() { return MPI_BYTE; }
	static int Length() { return sizeof(T); }
}; 
*/

template<class T> 
class MPITraits
{
public:
	static MPI_Datatype Type();
	static int Length() { return 1; }
}; 

template<class T> 
class MPITraits< std::complex<T> >
{
public:
	static MPI_Datatype Type() { return MPITraits<T>::Type(); }
	static int Length() { return 2; }
};

template<> inline MPI_Datatype MPITraits<unsigned char>::Type()
{
	return MPI_BYTE;
}

template<> inline MPI_Datatype MPITraits<int>::Type()
{
	return MPI_INTEGER;
}

template<> inline MPI_Datatype MPITraits<float>::Type()
{
	return MPI_FLOAT;
}

template<> inline MPI_Datatype MPITraits<double>::Type()
{
	return MPI_DOUBLE;
}

#endif


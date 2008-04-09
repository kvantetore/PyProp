
#include <mpi.h>

extern "C"
{
	void FORTRAN_NAME(bandedmatrixmultiply)(cplx* potential, int* potentialExtent0, double* scaling, cplx* source, int* sourceExtent0, cplx* dest, int* destExtent0, int* procCount0, int* procId0, int* communicator0);
}

void BandedMatrixMultiply_Wrapper(Array< cplx, 1 > potential, double scaling, Array< cplx, 1 > source, Array< cplx, 1 > dest)
{
	int potentialExtent0 = potential.extent(0);
	int sourceExtent0 = source.extent(0);
	int destExtent0 = dest.extent(0);

	int procId, procCount;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);
	MPI_Fint commHandle = MPI_Comm_c2f(MPI_COMM_WORLD);
	
	FORTRAN_NAME(bandedmatrixmultiply)(potential.data(), &potentialExtent0, &scaling, source.data(), &sourceExtent0, dest.data(), &destExtent0, &procCount, &procId, &commHandle);
}



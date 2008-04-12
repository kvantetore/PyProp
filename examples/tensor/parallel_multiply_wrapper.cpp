#include <mpi.h>

/*
extern "C"
{
	void FORTRAN_NAME(bandedmatrixmultiply)(cplx* potential, int* potentialExtent0, double* scaling, cplx* source, int* sourceExtent0, cplx* dest, int* destExtent0, int* globalSize0, int* bands0, int* procCount0, int* procId0, int* communicator0);
}

void BandedMatrixMultiply_Wrapper(Array< cplx, 1 > potential, double scaling, Array< cplx, 1 > source, Array< cplx, 1 > dest, int globalSize0, int bands0)
{
	int potentialExtent0 = potential.extent(0);
	int sourceExtent0 = source.extent(0);
	int destExtent0 = dest.extent(0);

	int procId, procCount;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);
	MPI_Fint commHandle = MPI_Comm_c2f(MPI_COMM_WORLD);
	
	FORTRAN_NAME(bandedmatrixmultiply)(potential.data(), &potentialExtent0, &scaling, source.data(), &sourceExtent0, dest.data(), &destExtent0, &globalSize0, &bands0, &procCount, &procId, &commHandle);
}
*/

extern "C"
{
void FORTRAN_NAME(multiplyoverlapmatrix)(cplx* overlap, int* overlapExtent0, int* overlapExtent1, cplx* source, int* sourceExtent0, int* sourceExtent1, int* sourceExtent2, cplx* dest, int* destExtent0, int* destExtent1, int* destExtent2);
}

void MultiplyOverlapMatrix_Wrapper(Array< cplx, 2 > overlap, Array< cplx, 3 > source, Array< cplx, 3 > dest)
{
	int overlapExtent0 = overlap.extent(0);
	int overlapExtent1 = overlap.extent(1);

	int sourceExtent0 = source.extent(0);
	int sourceExtent1 = source.extent(1);
	int sourceExtent2 = source.extent(2);

	int destExtent0 = dest.extent(0);
	int destExtent1 = dest.extent(1);
	int destExtent2 = dest.extent(2);

	FORTRAN_NAME(multiplyoverlapmatrix)(overlap.data(), &overlapExtent0, &overlapExtent1, source.data(), &sourceExtent0, &sourceExtent1, &sourceExtent2, dest.data(), &destExtent0, &destExtent1, &destExtent2);

}



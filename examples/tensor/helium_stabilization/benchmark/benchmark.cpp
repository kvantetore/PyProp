#include <complex>
#include <iostream>
#include <mpi.h>

#include "blitzhdf.h"

#define FORTRAN_NAME(x) x ## _
using namespace blitz;

//Generated code for wrapping fortran function
extern "C"
{
		void FORTRAN_NAME(tensorpotentialmultiply_simpd_simp_bandnh)(cplx* potential, int* potentialExtent0, int* potentialExtent1, int* potentialExtent2, double* scaling, cplx* source, int* sourceExtent0, int* sourceExtent1, int* sourceExtent2, cplx* dest, int* destExtent0, int* destExtent1, int* destExtent2, int* globalSize0, int* localMatrixIndex0, int* localMatrixIndex0Extent0, int* globalRow0, int* globalRow0Extent0, int* globalCol0, int* globalCol0Extent0, int* sendProc0, int* sendProc0Extent0, int* recvProcList0, int* recvProcList0Extent0, int* recvProcList0Extent1, int* recvLocalRowList0, int* recvLocalRowList0Extent0, int* recvLocalRowList0Extent1, int* recvCount0, int* recvCount0Extent0, cplx* recvTemp0, int* recvTemp0Extent0, int* recvTemp0Extent1, int* recvTemp0Extent2, cplx* sendTemp0, int* sendTemp0Extent0, int* sendTemp0Extent1, int* sendTemp0Extent2, int* pair1, int* pair1Extent0, int* pair1Extent1);
}
	
void TensorPotentialMultiply_SimpD_Simp_BandNH_Wrapper(Array< cplx, 3 > potential, double scaling, Array< cplx, 3 > source, Array< cplx, 3 > dest, int globalSize0, Array< int, 1 > localMatrixIndex0, Array< int, 1 > globalRow0, Array< int, 1 > globalCol0, Array< int, 1 > sendProc0, Array< int, 2 > recvProcList0, Array< int, 2 > recvLocalRowList0, Array< int, 1 > recvCount0, Array< cplx, 3 > recvTemp0, Array< cplx, 3 > sendTemp0, Array< int, 2 > pair1)
	{
		int potentialExtent0 = potential.extent(0);
		int potentialExtent1 = potential.extent(1);
		int potentialExtent2 = potential.extent(2);
		int sourceExtent0 = source.extent(0);
		int sourceExtent1 = source.extent(1);
		int sourceExtent2 = source.extent(2);
		int destExtent0 = dest.extent(0);
		int destExtent1 = dest.extent(1);
		int destExtent2 = dest.extent(2);
		int localMatrixIndex0Extent0 = localMatrixIndex0.extent(0);
		int globalRow0Extent0 = globalRow0.extent(0);
		int globalCol0Extent0 = globalCol0.extent(0);
		int sendProc0Extent0 = sendProc0.extent(0);
		int recvProcList0Extent0 = recvProcList0.extent(0);
		int recvProcList0Extent1 = recvProcList0.extent(1);
		int recvLocalRowList0Extent0 = recvLocalRowList0.extent(0);
		int recvLocalRowList0Extent1 = recvLocalRowList0.extent(1);
		int recvCount0Extent0 = recvCount0.extent(0);
		int recvTemp0Extent0 = recvTemp0.extent(0);
		int recvTemp0Extent1 = recvTemp0.extent(1);
		int recvTemp0Extent2 = recvTemp0.extent(2);
		int sendTemp0Extent0 = sendTemp0.extent(0);
		int sendTemp0Extent1 = sendTemp0.extent(1);
		int sendTemp0Extent2 = sendTemp0.extent(2);
		int pair1Extent0 = pair1.extent(0);
		int pair1Extent1 = pair1.extent(1);
		
		FORTRAN_NAME(tensorpotentialmultiply_simpd_simp_bandnh)(potential.data(), &potentialExtent0, &potentialExtent1, &potentialExtent2, &scaling, source.data(), &sourceExtent0, &sourceExtent1, &sourceExtent2, dest.data(), &destExtent0, &destExtent1, &destExtent2, &globalSize0, localMatrixIndex0.data(), &localMatrixIndex0Extent0, globalRow0.data(), &globalRow0Extent0, globalCol0.data(), &globalCol0Extent0, sendProc0.data(), &sendProc0Extent0, recvProcList0.data(), &recvProcList0Extent0, &recvProcList0Extent1, recvLocalRowList0.data(), &recvLocalRowList0Extent0, &recvLocalRowList0Extent1, recvCount0.data(), &recvCount0Extent0, recvTemp0.data(), &recvTemp0Extent0, &recvTemp0Extent1, &recvTemp0Extent2, sendTemp0.data(), &sendTemp0Extent0, &sendTemp0Extent1, &sendTemp0Extent2, pair1.data(), &pair1Extent0, &pair1Extent1);
}

using std::cout;
using std::endl;

std::string GetArgumentDatasetName(int argumentNumber)
{
	std::stringstream datasetName;
	datasetName << "/argument_" << argumentNumber;
	return datasetName.str();
}

template<class T, int Rank> 
blitz::Array<T, Rank> ReadArgumentArray(hid_t fileId, int argumentNumber)
{
	blitz::Array<T, Rank> array;
	hid_t datasetId = H5Dopen(fileId, GetArgumentDatasetName(argumentNumber).c_str());
	if (datasetId >= 0)
	{
		array.reference( ReadArray<T, Rank>(datasetId) );
		H5Dclose(datasetId);
	}
	return array;
}

template<class T> 
T ReadArgumentScalar(hid_t fileId, int argumentNumber)
{
	hid_t datasetId = H5Dopen(fileId, GetArgumentDatasetName(argumentNumber).c_str());
	T value = ReadScalar<T>(datasetId);
	H5Dclose(datasetId);
	return value;
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int procId, procCount;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);

	//open hdf file
	std::stringstream filename;
	filename << "potential_" << procId << ".h5";
	std::string filenameStr = filename.str();
	cout << "trying to open file " << filenameStr << endl;
	cout << "..." << endl;
	hid_t fileId = H5Fopen(filenameStr.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	//read arguments
	blitz::Array< cplx, 3 > potential = ReadArgumentArray< cplx, 3 >(fileId, 0	);
	double scaling = ReadArgumentScalar<double>(fileId, 1);
	blitz::Array< cplx, 3 > source = ReadArgumentArray< cplx, 3 >(fileId, 2	);
	blitz::Array< cplx, 3 > dest = ReadArgumentArray< cplx, 3 >(fileId, 3	);
	int globalSize0 = ReadArgumentScalar<int>(fileId, 4);
	blitz::Array< int, 1 > localMatrixIndex0 = ReadArgumentArray< int, 1 >(fileId, 5	);
	blitz::Array< int, 1 > globalRow0 = ReadArgumentArray< int, 1 >(fileId, 6	);
	blitz::Array< int, 1 > globalCol0 = ReadArgumentArray< int, 1 >(fileId, 7	);
	blitz::Array< int, 1 > sendProc0 = ReadArgumentArray< int, 1 >(fileId, 8	);
	blitz::Array< int, 2 > recvProcList0 = ReadArgumentArray< int, 2 >(fileId, 9	);
	blitz::Array< int, 2 > recvLocalRowList0 = ReadArgumentArray< int, 2 >(fileId, 10	);
	blitz::Array< int, 1 > recvCount0 = ReadArgumentArray< int, 1 >(fileId, 11	);
	blitz::Array< cplx, 3 > recvTemp0 = ReadArgumentArray< cplx, 3 >(fileId, 12	);
	blitz::Array< cplx, 3 > sendTemp0 = ReadArgumentArray< cplx, 3 >(fileId, 13	);
	blitz::Array< int, 2 > pair1 = ReadArgumentArray< int, 2 >(fileId, 14	);

	H5Fclose(fileId);

	MPI_Barrier(MPI_COMM_WORLD);
	cout << procId << " says hei" << endl;

	TensorPotentialMultiply_SimpD_Simp_BandNH_Wrapper(potential, scaling, source, dest, globalSize0, localMatrixIndex0, globalRow0, globalCol0, sendProc0, recvProcList0, recvLocalRowList0, recvCount0, recvTemp0, sendTemp0, pair1);

	cout << procId << " says på deg" << endl;


	MPI_Finalize();
}


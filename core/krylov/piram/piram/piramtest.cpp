#include "piram.h"

#include <iostream>
#include <fstream>
#include <strstream>
#include <mpi.h>

using namespace std;
using namespace blitz;
using namespace piram;

typedef TinyVector<int, 1> DataVector;
typedef Array<cplx, 1> DataArray;

const int MaxMatrixSize = 1000000;

struct MatrixElement
{
	int RowIndex;
	int ColIndex;

	double value;
};

int procId;
int procCount;
Array< MatrixElement, 1 > Matrix(MaxMatrixSize);
Array< cplx, 1 > GlobalIn;
Array< cplx, 1 > GlobalOut;
int matrixSize = 0;
int globalSize = 0;
int localSize = 0;

void LoadMatrix()
{
	ifstream file;
	file.open("qdotmatrix.dat", ios::in);

	int i=0;
	while (! file.eof() && file.is_open())
	{
		char line[256];
		file.getline(line, 256);
		//cout << "'" << line << "'" << endl;

		std::istrstream lineStream(line);
		double v0, v1;
		lineStream >> Matrix(i).RowIndex;
		lineStream >> Matrix(i).ColIndex;
		lineStream >> v0;
		lineStream >> v1;

		Matrix(i).value = v0;
		//cout << Matrix(i).RowIndex << " " << Matrix(i).ColIndex << " " << v0 << " " << v1 << endl;
		i++;
	}

	matrixSize = i-1;
	globalSize = Matrix(matrixSize-1).RowIndex+1;
	localSize = globalSize / procCount;
	if (globalSize % procCount != 0)
	{
		cout << "ERROR! globalSize (" << globalSize << ") is not div by procCount (" << procCount <<")" << endl;
	}

	GlobalIn.resize(globalSize);
	GlobalOut.resize(globalSize);
}


class OpFunctor : public OperatorFunctor<cplx>
{
public:
	typedef boost::shared_ptr< OpFunctor > Ptr;

	virtual void operator()(blitz::Array<cplx, 1> &localIn, blitz::Array<cplx, 1> &localOut) 
	{
		int localStart = localSize * procId;
		static int count = 0;
		
		//Allgather local in
		GlobalIn = 0;
		MPI_Allgather(localIn.data(), localSize*2, MPI_DOUBLE, GlobalIn.data(), localSize*2, MPI_DOUBLE, MPI_COMM_WORLD);
		//if (count % 10 == 0) cout << "Iteration " << count << endl;
		count++;
		//Perform matrix-vector multiply on global vector
		GlobalOut = 0;
		for (int i=0; i<matrixSize; i++)
		{
			int row = Matrix(i).RowIndex;
			int col = Matrix(i).ColIndex;
			double value = Matrix(i).value;
		
			GlobalOut(row) += value * GlobalIn(col);
			//GlobalOut(col) += value * GlobalIn(row);
		}
		
		//Assign the correct part to localOut
		localOut = GlobalOut( Range(localStart, localStart+localSize-1) );
	}
};

class SetupFunctor : public SetupResidualFunctor<cplx>
{
public:
	typedef boost::shared_ptr< SetupFunctor > Ptr;

	virtual void operator()(blitz::Array<cplx, 1> &residual) 
	{
		int procId, procCount;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		MPI_Comm_size(MPI_COMM_WORLD, &procCount);

		int MatrixSize = residual.extent(0);

		int localStart = procId * MatrixSize;
		cout << "proc " << procId <<": LocalStart = " << localStart << endl;

		ifstream file;
		file.open("init.dat", ios::in);

		int i=0;
		int localIndex = 0;
		while (! file.eof() && file.is_open())
		{
			
			char line[256];
			file.getline(line, 256);

			if (i>=localStart && i<localStart+MatrixSize)
			{
				std::istrstream lineStream(line);
				lineStream >> residual(localIndex);;
				localIndex++;
			}
			
			i++;
		}


		file.close();

	}

};




template class pIRAM<cplx>;
typedef pIRAM<cplx> pIRAMType;

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);

	LoadMatrix();

	OpFunctor::Ptr op = OpFunctor::Ptr( new OpFunctor() );
	SetupFunctor::Ptr setup= SetupFunctor::Ptr( new SetupFunctor() );


	pIRAMType piram;

	piram.MatrixOperator = op;
	piram.SetupResidual = setup;
	piram.UseRandomStart = true;

	piram.MaxRestartCount = 10;
	piram.MaxOrthogonalizationCount = 3;
	piram.EigenvalueCount = 5;
	piram.MatrixSize = localSize;
	piram.BasisSize = 50;
	piram.Tolerance = 0;

	piram.Setup();
	piram.Solve();
	piram.Postprocess();
}



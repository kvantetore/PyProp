#include <core/common.h>
#include <core/wavefunction.h>
#include <core/representation/representation.h>
#include <core/trilinos/pyprop_epetra.h>
#include <sstream>

template<int Rank>
class EpetraPotential
{
public:
	EpetraPotential() {}
	~EpetraPotential() {}

	typedef blitz::Array<cplx, Rank> ArrayType;

	void Setup(typename Wavefunction<Rank>::Ptr psi)
	{
		//Create Epetra map from wavefunction
		ProcessorMap = CreateWavefunctionEpetraMap<Rank>(psi);

		//Input and output vectors (views to wavefunctions)
		InVector = Epetra_Vector_Ptr( new Epetra_Vector(View, *ProcessorMap, (double*)psi->GetData().data()) );
		OutVector = Epetra_Vector_Ptr( new Epetra_Vector(View, *ProcessorMap, (double*)psi->GetData().data()) );
		
		//Get strides of the global Rank-dimensional array
		blitz::TinyVector<int, Rank> globalShape = psi->GetRepresentation()->GetFullShape();
		GlobalStrides(Rank-1) = 1;
		for (int rank=Rank-2; rank>=0; rank--)
		{
			GlobalStrides(rank) = GlobalStrides(rank+1)	* globalShape(rank+1);
		}

		//Convert tensor potential to Epetra FECrsMatrix
		Matrix = Epetra_FECrsMatrix_Ptr( new Epetra_FECrsMatrix(Copy, *ProcessorMap, 0) );

		//Matrix = CreateTensorPotentialEpetraMatrix<Rank>(*ProcessorMap, potentialData, basisPairs, globalStrides, cutoff, false);
	
		//Ignore warning from Epetra (every time Matrix is expanded a warning is issued - annoying)
		Matrix->SetTracebackMode(1);

		//Print some global matrix properties
	//	stringstream infoStr(stringstream::in | stringstream::out);
	//	infoStr << "\nEPETRA MATRIX INFO" << endl;
	//	infoStr << "  ||M||_F = " << Matrix->NormFrobenius() << endl;
	//	infoStr << "  Global nonzeros = " << Matrix->NumGlobalNonzeros() << endl;
	//	infoStr << "  Tensor potential size = " << potentialData.size() << endl;
	//	infoStr << endl;
	//
	//	cout << infoStr.str();
	}

	void AddTensorPotentialData(ArrayType potentialData, boost::python::list basisPairs, double cutOff)
	{
		CopyTensorPotentialToEpetraMatrix(Matrix, potentialData, basisPairs, GlobalStrides, cutOff);
	}

	void GlobalAssemble()
	{
		Matrix->GlobalAssemble();
	}

	void Multiply(typename Wavefunction<Rank>::Ptr srcPsi, typename Wavefunction<Rank>::Ptr destPsi)
	{
		//Update InVector to point to srcPsi data
		InVector->ResetView((double*)srcPsi->GetData().data());

		//Update OutVector to point to destPsi data
		OutVector->ResetView((double*)destPsi->GetData().data());

		//Multiply tensor potential on srcPsi
		Matrix->Multiply(false, *InVector, *OutVector);
	}


private:
	Epetra_Map_Ptr ProcessorMap;
	blitz::TinyVector<int, Rank> GlobalStrides;
	Epetra_FECrsMatrix_Ptr Matrix;
	Epetra_Vector_Ptr InVector;
	Epetra_Vector_Ptr OutVector;
};


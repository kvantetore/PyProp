#include <core/common.h>
#include <core/wavefunction.h>
#include <core/representation/representation.h>
#include <core/trilinos/pyprop_epetra.h>

#include <Ifpack.h>

typedef shared_ptr<Ifpack_Preconditioner> Ifpack_Preconditioner_Ptr;

template<int Rank>
class IfpackRadialPreconditioner
{
public:
	IfpackRadialPreconditioner() {}
	~IfpackRadialPreconditioner() {}

	typedef blitz::Array<cplx, Rank> ArrayType;

	void Setup(ArrayType vectorData, ArrayType potentialData, boost::python::list basisPairs, double cutoff)
	{
		Epetra_SerialComm comm;
		blitz::TinyVector<int, Rank> stride = vectorData.stride();
		Epetra_LocalMap map(vectorData.size()*2, 0, comm);
		Matrix = CreateTensorPotentialEpetraMatrix(map, potentialData, basisPairs, stride, cutoff);

		InVector = Epetra_Vector_Ptr( new Epetra_Vector(View, map, (double*)vectorData.data()) );
		OutVector = Epetra_Vector_Ptr( new Epetra_Vector(map) );
		
		// create the preconditioner. For valid PrecType values,
		// please check the documentation
		string PrecType = "ILU"; // incomplete LU
		int OverlapLevel = 0; // must be >= 0. If Comm.NumProc() == 1,
		                      // it is ignored.
		
		Ifpack Factory;
		Preconditioner = Ifpack_Preconditioner_Ptr( Factory.Create(PrecType, Matrix.get(), OverlapLevel) );
		assert(Preconditioner != 0);
		
		// specify parameters for ILU
		Parameters.set("fact: drop tolerance", 1e-9);
		Parameters.set("fact: level-of-fill", 1);
		// the combine mode is on the following:
		// "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
		// Their meaning is as defined in file Epetra_CombineMode.h   
		//Parameters.set("schwarz: combine mode", "Add");
		// sets the parameters
		Preconditioner->SetParameters(Parameters);
		
		// initialize the preconditioner. At this point the matrix must
		// have been FillComplete()'d, but actual values are ignored.
		Preconditioner->Initialize();
		
		// Builds the preconditioners, by looking for the values of 
		// the matrix.
		Preconditioner->Compute();
	}

	void Solve(ArrayType data)
	{
		//Update InVector to point to the RHS data
		InVector->ResetView((double*)data.data());

		//Solve preconditioner
		Preconditioner->ApplyInverse(*InVector, *OutVector);

		//Copy data from out vector back to array
		OutVector->ExtractCopy((double*)data.data());
	}

private:
	Epetra_FECrsMatrix_Ptr Matrix;
	Epetra_Vector_Ptr InVector;
	Epetra_Vector_Ptr OutVector;
	Ifpack_Preconditioner_Ptr Preconditioner;
	Teuchos::ParameterList Parameters;
};



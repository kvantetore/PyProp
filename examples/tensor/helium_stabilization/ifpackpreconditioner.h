#include <core/common.h>
#include <core/wavefunction.h>
#include <core/representation/representation.h>
#include <core/trilinos/pyprop_epetra.h>
#include <core/trilinos/pyprop_tpetra.h>

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
		//string PrecType = "Amesos"; // incomplete LU
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

#include <AnasaziConfigDefs.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBlockDavidsonSolMgr.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziLOBPCGSolMgr.hpp>
#include <AnasaziBasicOutputManager.hpp>
#include <AnasaziEpetraAdapter.hpp>
#include <AnasaziTpetraAdapter.hpp>
#include <Epetra_InvOperator.h>
#include <Teuchos_ArrayRCP.hpp>


template<int Rank>
class PypropOperator: public Tpetra_Operator
{
public:
	typedef shared_ptr< PypropOperator<Rank> > Ptr;
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;


	PypropOperator(object callback, typename Wavefunction<Rank>::Ptr srcPsi, typename Wavefunction<Rank>::Ptr dstPsi) : 
		Callback(callback),
		SourcePsi(srcPsi),
		DestPsi(dstPsi)
	{
		WavefunctionMap = CreateWavefunctionTpetraMap<Rank>(srcPsi);
	}

	//---------- Tpetra Operator Interface -----------------
	
    virtual int SetUseTranspose(bool UseTranspose)
	{
		return -1;
	}
  
    virtual void apply(const Tpetra_MultiVector& Xconst, Tpetra_MultiVector& Y, Teuchos::ETransp mode = Teuchos::NO_TRANS) const 
	{
		Tpetra_MultiVector &X = const_cast<Tpetra_MultiVector&>(Xconst);
		
		//remember original data buffers
		DataArray oldData = SourcePsi->GetData();
		DataArray oldTempData = DestPsi->GetData();

		//get shapes
		DataVector shape = SourcePsi->GetData().shape();
		DataVector stride = SourcePsi->GetData().stride();
	
		//Get local data from epetra vectors
		int ldaIn = 0;
		int ldaOut = 0;
		cplx *in;
		cplx *out;
		X.extractView1D(in, ldaIn);
		Y.extractView1D(out, ldaOut);

		//cout << "Input = " << X << endl;

		//Iterate over each vector
		for (int i=0; i<Y.numVectors(); i++)
		{
			//put epetra data into blitz arrays
			DataArray inData((cplx*)in, shape, stride, blitz::neverDeleteData);
			DataArray outData((cplx*)out, shape, stride, blitz::neverDeleteData);
			outData = 0;
			
			//Set psi and tempPsi to point to correct vectorsbuffers
			SourcePsi->SetData(inData);
			DestPsi->SetData(outData);
			
			Callback(SourcePsi, DestPsi);

			in += ldaIn;
			out += ldaOut;
		}
	
		//cout << "Output = " << Y << endl;
		
		//Restore the former buffers
		SourcePsi->SetData(oldData);
		DestPsi->SetData(oldTempData);
	}

    //! Returns the Tpetra_Map object associated with the domain of this operator.
    virtual const Tpetra_Map & getDomainMap() const 
	{
		return *WavefunctionMap;
	}

    //! Returns the Tpetra_Map object associated with the range of this operator.
    virtual const Tpetra_Map & getRangeMap() const 
	{
		return *WavefunctionMap;
	}


private:
	object Callback;
	typename Wavefunction<Rank>::Ptr SourcePsi;
	typename Wavefunction<Rank>::Ptr DestPsi;

	Tpetra_Map_Ptr WavefunctionMap;
	Tpetra_Comm_Ptr TpetraComm;
};

using Teuchos::RCP;
using Teuchos::rcp;

template<int Rank>
class AnasaziSolver 
{
public:
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;

	typedef cplx SC;
	typedef Tpetra_MultiVector MV;
	typedef Tpetra_Operator OP;
	typedef Anasazi::MultiVecTraits<SC, MV> MVT;
	
	typedef Anasazi::BasicEigenproblem<SC, MV, OP> EigenProb;
	typedef RCP<EigenProb> EigenProb_Ptr;
	typedef Anasazi::SolverManager<SC, MV, OP> SolverMgr;
	typedef Anasazi::Eigensolution<SC, MV> EigenSol;

private:
	//Temporary member variables
	typename Wavefunction<Rank>::Ptr Psi;
	typename Wavefunction<Rank>::Ptr TempPsi;
	EigenProb_Ptr Problem;
	

	int EigenvalueCount;
	int BlockSize;
	int BlockCount;
	int MaxIterationCount; 
	double Tolerance;
	std::string Which;
	bool Locking;
	double LockingTolerance;
	std::string Method;

	bool GeneralizedEigenvalueProblem;
	std::string OrthoMethod;

	bool PrintWarnings;               /*!< Internal warnings */
	bool PrintIterationDetails;       /*!< Approximate eigenvalues, errors */
	bool PrintOrthoDetails;           /*!< Orthogonalization/orthonormalization details */
	bool PrintFinalSummary;           /*!< Final computational summary */
	bool PrintTimingDetails;          /*!< Timing details */
	bool PrintStatusTestDetails;      /*!< Status test details */
	bool PrintDebug;                  /*!< Debugging information */
	bool PrintErrors;                 /*!< Errors [ always printed ] */


public:
	AnasaziSolver() {}
	virtual ~AnasaziSolver() {}

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("krylov_eigenvalue_count", EigenvalueCount);
		config.Get("krylov_block_size", BlockSize);
		config.Get("krylov_block_count", BlockCount);
		config.Get("krylov_max_iteration_count", MaxIterationCount);
		config.Get("krylov_tolerance", Tolerance);
		config.Get("krylov_which", Which);
		config.Get("krylov_locking", Locking);
		config.Get("krylov_locking_tolerance", LockingTolerance);
		config.Get("krylov_method", Method);
		config.Get("generalized_eigenvalue_problem", GeneralizedEigenvalueProblem);
		config.Get("orthogonalization", OrthoMethod);
		config.Get("print_warnings" ,PrintWarnings );
		config.Get("print_iteration_details" ,PrintIterationDetails );
		config.Get("print_ortho_details" ,PrintOrthoDetails );
		config.Get("print_final_summary" ,PrintFinalSummary );
		config.Get("print_timing_details" ,PrintTimingDetails );
		config.Get("print_statustest_details" ,PrintStatusTestDetails );
		config.Get("print_debug" ,PrintDebug );
		config.Get("print_errors" ,PrintErrors );
	}

	void Setup(const typename Wavefunction<Rank>::Ptr psi)
	{
	}

	void Solve(object operatorCallback, object preconditionerCallback, object overlapCallback, typename Wavefunction<Rank>::Ptr psi, typename Wavefunction<Rank>::Ptr tempPsi)
	{
		//Create EpetraOperator of callback
		RCP< PypropOperator<Rank> > op = Teuchos::rcp( new PypropOperator<Rank>(operatorCallback, psi, tempPsi) );
		RCP< PypropOperator<Rank> > overOp = Teuchos::rcp( new PypropOperator<Rank>(overlapCallback, psi, tempPsi) );
		RCP< PypropOperator<Rank> > precOp = Teuchos::rcp( new PypropOperator<Rank>(preconditionerCallback, psi, tempPsi) );

		//Create Map from wavefunction
		Tpetra_Map_Ptr map = CreateWavefunctionTpetraMap<Rank>(psi);

		//Create input vector
		RCP<MV> inputVector = rcp( new MV(*map, BlockSize) );

		//Create problem
		Problem = rcp( new EigenProb(op, inputVector) ); 
		Problem->setHermitian(true);
		if (GeneralizedEigenvalueProblem)
		{
			Problem->setM(overOp);
		}
		if (preconditionerCallback)
		{
			Problem->setPrec(precOp);
		}
		Problem->setNEV(EigenvalueCount);
		if (!Problem->setProblem())
		{
			throw std::runtime_error("Error in setProblem()");
		}

		//Set parameters for solver
		Teuchos::ParameterList params;
		params.set("Which", Which);
		params.set("Block Size", BlockSize);
		params.set("Num Blocks", BlockCount);
		params.set("Maximum Restarts", MaxIterationCount);
		params.set("Convergence Tolerance", Tolerance);
		params.set("Orthogonalization", "DGKS");
		params.set("Verbosity", 
				(PrintWarnings ? Anasazi::Warnings & PrintWarnings : 0) |
				(PrintIterationDetails ? Anasazi::IterationDetails : 0 ) |
				(PrintOrthoDetails ? Anasazi::OrthoDetails : 0 ) |
				(PrintFinalSummary ? Anasazi::FinalSummary : 0 ) |
				(PrintTimingDetails ? Anasazi::TimingDetails : 0 ) |
				(PrintStatusTestDetails ? Anasazi::StatusTestDetails : 0 ) |
				(PrintDebug ? Anasazi::Debug : 0 ) |
				(PrintErrors ? Anasazi::Errors : 0 )
			);
		params.set("Use Locking", Locking);
		params.set("Locking Tolerance", LockingTolerance);

		RCP< SolverMgr > manager;
		if (Method == std::string("Davidson"))
		{
			//Create solver manager
			manager = rcp( new Anasazi::BlockDavidsonSolMgr<SC, MV, OP>(Problem, params) );
		}
		else if (Method == std::string("KrylovSchur"))
		{
			params.set("Orthogonalization", OrthoMethod);
			manager = rcp( new Anasazi::BlockKrylovSchurSolMgr<SC, MV, OP>(Problem, params) );
		}
		else if (Method == std::string("LOBPCG"))
		{
			params.set("Full Ortho", true);
			manager = rcp( new Anasazi::LOBPCGSolMgr<SC, MV, OP>(Problem, params) );
		}
		else
		{
			throw std::runtime_error("Unknown krylov_method (check case)");
		}
		manager->solve();

		//Get the solution
		EigenSol solution = Problem->getSolution();

		//print eigenvalues
		if (psi->GetRepresentation()->GetDistributedModel()->ProcId == 0)
		{
			cout << "Converged Eigenvalues = " << solution.numVecs << endl;
			for (int i=0; i<solution.numVecs; i++) 
			{
				cout
					<<	std::setw(16) << solution.Evals[i].realpart
					<<	std::endl;
			}
		}
	}

	void virtual SetupResidual(blitz::Array<cplx, 1> &residual)
	{
		throw std::runtime_error("SetupResidual not implemented");
	}

	/*
	 * Returns the converged eigenvalues
	 */
	blitz::Array<cplx, 1> GetEigenvalues()
	{
		EigenSol solution = Problem->getSolution();
		int numVec = solution.numVecs;
		blitz::Array<cplx, 1> ev(numVec);
		for	(int i=0; i<numVec; i++)
		{
			ev(i) = cplx(solution.Evals[i].realpart, solution.Evals[i].imagpart);
		}

		return ev;
	}

	/*
	 * Returns the _converged_ eigenvectors, a N by M matrix, where
	 * N == GetEigenvalueCount(), and M is the size of the wavefunction
	 */
	blitz::Array<cplx, 1> GetEigenvector(int eigenvectorIndex)
	{
		int lda;
		cplx* data;
		EigenSol solution = Problem->getSolution();
		solution.Evecs->extractView1D(data, lda);
		data += lda*eigenvectorIndex;
		blitz::TinyVector<int, 1> stride = 1;
		blitz::TinyVector<int, 1> shape = solution.Evecs->myLength();
		return blitz::Array<cplx, 1>(data, shape, stride, blitz::neverDeleteData);
	}

	double EstimateMemoryUsage(int matrixSize, int basisSize)
	{
		return -1.;
	}

	blitz::Array<double, 1> GetErrorEstimates()
	{
		throw std::runtime_error("Cannot get error estimates yet");
	}


	blitz::Array<double, 1> GetConvergenceEstimates()
	{
		throw std::runtime_error("Cannot get error estimates yet");
	}


	/*
	 * Returns the number of converged eigenvalues
	 * <= EigenvalueCount
	 */
	int GetEigenvalueCount()
	{
		return Problem->getSolution().numVecs;
	}
	
	/*
	 * Returns the number of restarting steps performed
	 * <= MaxIterationCount
	 */
	int GetRestartCount()
	{
		throw std::runtime_error("Cannot get restart count yet");
	}

	/* 
	 * Returns the nuber of matrix-vector operations performed
	 */
	int GetOperatorCount()
	{
		throw std::runtime_error("Cannot get operator count yet");
	}

	/*
	 * Returns the number of re-orthogonalization steps performed
	 */
	int GetOrthogonalizationCount()
	{
		throw std::runtime_error("Cannot get orthogonalization count yet");
	}

};





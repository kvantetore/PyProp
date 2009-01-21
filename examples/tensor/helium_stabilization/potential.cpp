#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

using namespace blitz;

/*
 * Radial Kinetic Energy
 */
template<int Rank>
class KineticEnergyPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double mass;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("mass", mass);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		return - 1. / (2. * mass);
	}

};

/* 
 * Coulomb Potential
 */
template<int Rank>
class CoupledSphericalCoulombPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	double Z;
	int RadialRank1;
	int RadialRank2;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("z", Z);
		config.Get("radial_rank1", RadialRank1);
		config.Get("radial_rank2", RadialRank2);
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r1 = pos(RadialRank1);
		double r2 = pos(RadialRank2);
		return - Z/r1 - Z/r2;
	}
};


template<int Rank>
class ComplexAbsorbingPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	int radialRank1;
	int radialRank2;
	double scalingReal;
	double scalingImag;
	double factorReal;
	double factorImag;
	double absorberStart;
	double absorberLength;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("radial_rank_1", radialRank1);
		config.Get("radial_rank_1", radialRank2);
		config.Get("absorber_start", absorberStart);
		config.Get("absorber_length", absorberLength);
		config.Get("scaling_real", scalingReal);
		config.Get("scaling_imag", scalingImag);
		config.Get("factor_real", factorReal);
		config.Get("factor_imag", factorImag);
	}

	/*
	 * Called for every grid point 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r1 = pos(radialRank1);
		double r2 = pos(radialRank2);
		double r = std::sqrt(r1*r1 + r2*r2);
		cplx V = 0;
		if (r > absorberStart)
		{
			double curLength = (r - absorberStart) / absorberLength;
			double Vr = factorReal * std::pow(curLength, scalingReal);
			double Vi = factorImag * std::pow(curLength, scalingImag);
			V = cplx(Vr , Vi);
		}
		return V;
	}
};


template<int Rank>
class OverlapPotential : public PotentialBase<Rank>
{
public:
        //Required by DynamicPotentialEvaluator
        cplx TimeStep;
        double CurTime;


        void ApplyConfigSection(const ConfigSection &config)
        {
        }

        inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
        {
                return 1.0;
        }
};


template<int Rank>
class BoxNormPotential : public PotentialBase<Rank>
{
public:
		int RadialRank1;
		int RadialRank2;
		int BoxSize;

        void ApplyConfigSection(const ConfigSection &config)
        {
			config.Get("radial_rank_1", RadialRank1);
			config.Get("radial_rank_2", RadialRank2);
			config.Get("box_size", BoxSize);
        }

        inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
        {
			double r1 = pos(RadialRank1);
			double r2 = pos(RadialRank2);
			double r = std::sqrt(r1 * r1 + r2 * r2);

			if (r < BoxSize)
				return 1.0;
			else 
				return 0.0;
        }
};

#include <slu_zdefs.h>
template<int Rank>
class SuperLUSolver
{
public:
	typedef blitz::Array<cplx, 1> VectorType;
	typedef blitz::Array<int, 1> VectorTypeInt;
	typedef blitz::Array<double, 1> VectorTypeReal;

	bool Initialized;
	VectorType MatrixData;
	blitz::Array<cplx, Rank> TempData;
	VectorTypeInt RowIndices;
	VectorTypeInt ColStartIndices;
	VectorTypeInt PermutationRow;
	VectorTypeInt PermutationCol;
	
	SuperMatrix A;
	SuperMatrix L; 
	SuperMatrix U;
	VectorTypeReal R;
	VectorTypeReal C;
	VectorTypeInt ETree;
	superlu_options_t Options;
    SuperLUStat_t Statistics;
	char Equed[2];
	

	SuperLUSolver()
	{
		Initialized = false;
	}

	~SuperLUSolver()
	{
		if (Initialized)
		{
			SUPERLU_FREE(A.Store);
			SUPERLU_FREE(L.Store);
			SUPERLU_FREE(U.Store);
		}
	}

	void Setup(int matrixSize, VectorType matrixData, VectorTypeInt rowIndices, VectorTypeInt colStartIndices)
	{
		/* Set the default input options:
		options.Fact = DOFACT;
		options.Equil = YES;
		options.ColPerm = COLAMD;
		options.DiagPivotThresh = 1.0;
		options.Trans = NOTRANS;
		options.IterRefine = NOREFINE;
		options.SymmetricMode = NO;
		options.PivotGrowth = NO;
		options.ConditionNumber = NO;
		options.PrintStat = YES;
		*/
		//Options.ColPerm = NATURAL;
		//Options.DiagPivotThresh = 0.1;
		//Options.SymmetricMode = YES;

		cout << "Starting SuperLU Setup" << endl;
		set_default_options(&Options);

		//Copy matrix data so it won't go out of scope
		MatrixData.reference( matrixData.copy() );
		RowIndices.reference( rowIndices.copy() );
		ColStartIndices.reference( colStartIndices.copy() );
		//We need permutations of rows and cols
		PermutationRow.resize( matrixSize*2 );
		PermutationCol.resize( matrixSize*2 );
		R.resize( matrixSize );
		C.resize( matrixSize );
		ETree.resize( matrixSize );

		int N = matrixSize;
		int nnz = MatrixData.size();

		//Create compressed col matrix from matrix data
		doublecomplex* data = (doublecomplex*) MatrixData.data();
    	zCreate_CompCol_Matrix(&A, N, N, nnz, data, RowIndices.data(), ColStartIndices.data(), SLU_NC, SLU_Z, SLU_GE);

   		//Create an empty right hand side to create factorization without solving anything
		int nrhs = 0;
		doublecomplex* rhs = 0;
		SuperMatrix B, X;
    	zCreate_Dense_Matrix(&B, N, nrhs, rhs, N, SLU_DN, SLU_Z, SLU_GE);
    	zCreate_Dense_Matrix(&X, N, nrhs, rhs, N, SLU_DN, SLU_Z, SLU_GE);

		//Factorize A
		int lwork = 0;
		doublecomplex *work = 0;

		double rpg; //reciprocal growth
		double ferr; //err[nrhs]
		double berr; //err[nrhs]
		double rcond;
		mem_usage_t memUsage;

		int info;
    	StatInit(&Statistics);
		Options.Fact = DOFACT;
		zgssvx(&Options, &A, PermutationCol.data(), PermutationRow.data(), ETree.data(), Equed, 
			R.data(), C.data(), &L, &U, work, lwork, &B, &X, &rpg, &rcond, &ferr, &berr,
           &memUsage, &Statistics, &info);

		if (info != 0)
		{
			cout << "Error from SuperLU: " << info << endl;
			throw std::runtime_error("Could not factorize matrix with superlu");
		}

		SCformat* Lstore = (SCformat *) L.Store;
		NCformat* Ustore = (NCformat *) U.Store;
		cout << "SuperLU-factorization of matrix returned " << info << endl;
		cout << "    Nonzero elemts in L = " << Lstore->nnz << endl;
		cout << "    Nonzero elemts in U = " << Ustore->nnz << endl;
		cout << "    Fill factor         = " << (double)(Lstore->nnz + Ustore->nnz - N) / (double) nnz << endl;

    	//Destroy temporary matrix B
		SUPERLU_FREE(B.Store);
		SUPERLU_FREE(X.Store);

		cout << "Finishing SuperLU Setup" << endl;

		Initialized = true;
	}

	void Solve(blitz::Array<cplx, Rank> rhs)
	{
		//allocate Temp-array if we don't already have it
		if (rhs.size() != TempData.size())
		{
			TempData.resize(rhs.shape());
		}

		//Setup right hand side in correct format
   		int nrhs = 1;
		doublecomplex* rhsData = (doublecomplex*)rhs.data();
		doublecomplex* tempData = (doublecomplex*)TempData.data();
		SuperMatrix B, X;
	 	zCreate_Dense_Matrix(&B, rhs.size(), nrhs, rhsData, rhs.size(), SLU_DN, SLU_Z, SLU_GE);
	 	zCreate_Dense_Matrix(&X, rhs.size(), nrhs, tempData, rhs.size(), SLU_DN, SLU_Z, SLU_GE);

		//Solve
		int lwork = 0;
		doublecomplex *work = 0;
		double rpg; //reciprocal growth
		double ferr; //err[nrhs]
		double berr; //err[nrhs]
		double rcond;
		mem_usage_t memUsage;

		int info;
		Options.Fact = FACTORED;
		zgssvx(&Options, &A, PermutationCol.data(), PermutationRow.data(), ETree.data(), Equed, 
			R.data(), C.data(), &L, &U, work, lwork, &B, &X, &rpg, &rcond, &ferr, &berr,
           &memUsage, &Statistics, &info);


		if (info != 0)
		{
			cout << "Error from SuperLU: " << info << endl;
			throw std::runtime_error("Could not factorize matrix with superlu");
		}

		//Copy data from temp data storage to overwrite input data
		rhs = TempData;

		//free temp datastructures
		SUPERLU_FREE(B.Store);
	}
};

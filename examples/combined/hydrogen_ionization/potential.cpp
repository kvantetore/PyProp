#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class LaserPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double PulseDuration;
	double Frequency;
	double Amplitude;

	//Calculated parameters
	double convolutionFrequency;
	double currentAmplitude;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("pulse_duration", PulseDuration);
		config.Get("frequency", Frequency);
		config.Get("amplitude", Amplitude);

		convolutionFrequency = M_PI / PulseDuration;
	}

	/*
	 * Called once every timestep
	 */
	void CurTimeUpdated()
	{
		if (CurTime > PulseDuration)
		{
			currentAmplitude = 0;
		}
		else
		{
			currentAmplitude = Amplitude;
			currentAmplitude *= sqr(sin(CurTime * convolutionFrequency));
			currentAmplitude *= cos(CurTime * Frequency);
		}
	}

	/*
	 * Called for every grid point at every time step. 
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = std::abs(pos(0));
		double theta = pos(1);

		double z = r * cos(theta);
		return currentAmplitude * z;
	}
};

template<int Rank>
class CoulombPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double charge;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("charge", charge);
	}

	/*
	 * Called for every grid point at every time step. 
	 *
	 * Some general tips for max efficiency:
	 * - If possible, move static computations to ApplyConfigSection.
	 * - Minimize the number of branches ("if"-statements are bad)
	 * - Minimize the number of function calls (sin, cos, exp, are bad)
	 * - Long statements can confuse the compiler, consider making more 
	 *   simpler statements
	 */
	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double r = std::abs(pos(0));
		if (r < 10e-5)
		{
			return 0;
		}
		double V = charge / r;

		return V;
	}
};

#include <core/finitediff/cranknicholsonpropagator.h>
#include <core/representation/combinedrepresentation.h>

template <int Rank>
class FiniteDifferenceSolver
{
public:
	typedef blitz::Array<cplx, 2> MatrixType;
	typedef blitz::Array<int, 1> IntVectorType;
	typedef blitz::Array<cplx, Rank> PotentialType;

	blitz::Array<MatrixType, 1> MatrixData;
	PotentialType PotentialData;
	blitz::Array<IntVectorType, 1> PivotData;
	int DifferenceRank;

	FiniteDifferenceSolver() {}
	~FiniteDifferenceSolver() {}

	void AddPotential(PotentialType potential)
	{
		if (PotentialData.size() == 0)
		{
			PotentialData.resize(potential.shape());
		}

		PotentialData += potential;
	}

	void Setup(typename Wavefunction<Rank>::Ptr psi, int differenceRank, int differenceOrder, cplx scalingI, cplx scalingH)
	{
		DifferenceRank = differenceRank;

		if (PotentialData.size() == 0)
		{
			throw std::runtime_error("No potentials added");
		}

		//Set up big matrix with correct number of matrices: 
		//  #elements before differenceRank * #elements after differenceRank
		int numMatrices = PotentialData.size() / PotentialData.extent(differenceRank);
		MatrixData.resize(numMatrices);
		PivotData.resize(numMatrices);

		//Get representation
		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast<CmbRepr>(psi->GetRepresentation());

		//Only orthogonal basises for now
		if (!repr->IsOrthogonalBasis(DifferenceRank))
		{
			throw std::runtime_error("DifferenceRank must be an orthogonal basis");
		}

		//Get Difference Matrix
		blitz::Array<double, 1> globalGrid = repr->GetGlobalGrid(differenceRank);
		FiniteDifferenceHelper fd;
		fd.Setup(globalGrid, differenceOrder);
		MatrixType laplacianBlasBanded = fd.SetupLaplacianBlasBanded();
		MatrixType laplacianLapackBanded = ConvertMatrixBlasBandedToLapackBanded(laplacianBlasBanded);

		double mass = 1.0;

		//Setup a matrix for all other grid points 
		blitz::Array<cplx, 3> potential3d = MapToRank3(PotentialData, differenceRank, 1);
		int matrixIndex = 0;
		for (int outerIndex=0; outerIndex<potential3d.extent(0); outerIndex++)
		{
			for (int innerIndex=0; innerIndex<potential3d.extent(2); innerIndex++)
			{
				//Create a view of the current matrix
				PivotData(matrixIndex).resize(laplacianLapackBanded.extent(0));
				MatrixData(matrixIndex).resize(laplacianLapackBanded.shape());
				MatrixType matrix = MatrixData(matrixIndex);
				IntVectorType pivots = PivotData(matrixIndex);
				blitz::Array<cplx, 1> diagonal = GetDiagonalViewLapackBanded(matrix);
				//Create a view of the current potential
				blitz::Array<cplx, 1> potentialSlice = potential3d(outerIndex, blitz::Range::all(), innerIndex);

				//Set matrix = (scalingI * I + scalingH + H) 
				//where I is the identity matrix and H is the laplacian plus potentials
				//Set laplacian
				matrix = laplacianLapackBanded;
				matrix *= -0.5 / mass;
				//add potentials
				diagonal += potentialSlice;
				//scale with the scaling coefficient
				matrix *= scalingH;
				//add identity * scaling
				diagonal += scalingI;

				//Create LU factorization
			 	lapack.CalculateLUFactorizationBanded(matrix, pivots);	
				
				matrixIndex++;
			}
		}
	}

	void Solve(typename Wavefunction<Rank>::Ptr psi)
	{
		//Get representation
		typedef CombinedRepresentation<Rank> CmbRepr;
		typename CmbRepr::Ptr repr = boost::static_pointer_cast<CmbRepr>(psi->GetRepresentation());

		//Only orthogonal basises for now
		blitz::Array<cplx, 3> psi3d = MapToRank3(psi->GetData(), DifferenceRank, 1);
		int rankSize = psi->GetData().extent(DifferenceRank);
		if (temp.extent(0) != rankSize)
		{
			temp.resize(rankSize);
		}

		//Solve
		int matrixIndex = 0;
		for (int outerIndex=0; outerIndex<psi3d.extent(0); outerIndex++)
		{
			for (int innerIndex=0; innerIndex<psi3d.extent(2); innerIndex++)
			{
				//Create a view of the current matrix
				MatrixType matrix = MatrixData(matrixIndex);
				IntVectorType pivots = PivotData(matrixIndex);
				
				//Create a copy of the current psi to make it unit strided
				temp = psi3d(outerIndex, blitz::Range::all(), innerIndex);

				//Solve for psi
				lapack.SolveBandedFactored(matrix, pivots, temp);

				//copy solution back into psi
				psi3d(outerIndex, blitz::Range::all(), innerIndex) = temp;

				matrixIndex++;
			}
		}
	}


private:
		typedef blitz::linalg::LAPACK<cplx> LAPACK;
		blitz::Array<cplx, 1> temp;
		LAPACK lapack;
};
	


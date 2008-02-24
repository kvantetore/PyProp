#include "bsplinepropagator.h"
#include "../representation/representation.h"
#include "../../utility/blitzblas.h"
#include "../../utility/blitztricks.h"
#include "../../utility/blitzlapack.h"
#include <cmath>


namespace BSpline
{

template<int Rank>
void Propagator<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	config.Get("mass", Mass);
	cout << "BSplinePropagator: Mass = " << Mass << endl;
}


template<int Rank>
void Propagator<Rank>::Setup(const cplx &dt, const Wavefunction<Rank> &psi, BSpline::Ptr bsplineObject, int rank)
{
	using namespace blitz;

	//Get b-spline object from wavefunction
	BSplineObject = bsplineObject;

	//Set class parameters
	PropagateRank = rank;

	int N = psi.GetRepresentation()->GetFullShape()(rank);

	//create some temporary arrays
	Array<cplx, 2> HamiltonianMatrixSetup;

	//Call setup routine to calculate Hamiltonian matrix
	SetupHamiltonianMatrix(HamiltonianMatrixSetup);

	//Obtain eigenvectors for HamiltonianMatrixSetup
	ComputeHamiltonianEigenvectors(HamiltonianMatrixSetup);

	//Create propagation matrix
	PropagationMatrix.resize(Eigenvectors.shape());

	/*
	 * Compute exponential of hamilton matrix H (-> propagation matrix)
	 * This is done using right and left eigenvectors of S**-1 * H,
	 * (X and Y**H respectively). Since X**-1 = Y**H, we have:
	 *
	 *     exp(H) = X * exp(L) * Y**H
	 *
	 * Eigenvector matrices are stored in col-major (FORTRAN) order.
	 */
	cplx expDot = 0;
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			expDot = 0;
			for (int k=0; k<N; k++)
			{
				expDot += Eigenvectors(k,i) * exp(-I * dt * Eigenvalues(k)) * conj(EigenvectorsInverse(k,j));
				//expDot += Eigenvectors(k,i) * conj(Eigenvectors(k,j));
			}
			PropagationMatrix(i,j) = expDot;
		}
	}

	//Allocate temp data
	TempData.resize(Eigenvalues.extent(0));
}

template<int Rank>
void Propagator<Rank>::AdvanceStep(Wavefunction<Rank> &psi)
{
	using namespace blitz;

	//Map the data to a 3D array, where the radial part is the 
	//middle rank
	Array<cplx, 3> data3d = MapToRank3(psi.Data, PropagateRank, 1);

	//Propagate the 3D array
	ApplyMatrix(PropagationMatrix, data3d);
}


template<int Rank>
void Propagator<Rank>::MultiplyHamiltonian(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi)
{
	using namespace blitz;

	//Map the data to a 3D array, where the radial part is the 
	//middle rank
	Array<cplx, 3> srcData = MapToRank3(srcPsi.Data, PropagateRank, 1);
	Array<cplx, 3> dstData = MapToRank3(dstPsi.Data, PropagateRank, 1);
	Array<cplx, 2> matrix = GetHamiltonianMatrix();

	//iterate over the array direction which is not propagated by this
	//propagator
	Array<cplx, 1> temp(srcData.extent(1));
	for (int i=0; i<srcData.extent(0); i++)
	{
		for (int j=0; j<srcData.extent(2); j++)
		{
			Array<cplx, 1> v = srcData(i, Range::all(), j);
			Array<cplx, 1> w = dstData(i, Range::all(), j);

			MatrixVectorMultiply(matrix, v, temp);
			w += temp;
		}
	}
}


template<int Rank>
void Propagator<Rank>::ApplyMatrix(const blitz::Array<cplx, 2> &matrix, blitz::Array<cplx, 3> &data)
{
	using namespace blitz;

	Array<cplx, 1> temp = TempData;
			
	//iterate over the array direction which is not propagated by this
	//propagator
	for (int i=0; i<data.extent(0); i++)
	{
		for (int j=0; j<data.extent(2); j++)
		{
			/* v is a view of a slice of the wave function
			 * temp is a temporary copy of that slice
			 */
			Array<cplx, 1> v = data(i, Range::all(), j);
			temp = v; 
			MatrixVectorMultiply(matrix, temp, v);
		}
	}
}


/*
 * Following are some internal functions to set up and diagonalize
 * the hamiltonian matrix in the b-spline basis.
 */
template<int Rank>
void Propagator<Rank>::SetupHamiltonianMatrix(blitz::Array<cplx, 2> &HamiltonianMatrix)
{
	int k = BSplineObject->MaxSplineOrder;
	int N = BSplineObject->NumberOfBSplines;
	HamiltonianMatrix.resize(N, N);
	HamiltonianMatrix = 0;

	OverlapMatrix.resize(N,N);
	OverlapMatrix = 0;

	for (int i = 0; i < N; i++)
	{
		int jMax = std::min(i + k, N);
		for (int j = i; j < jMax; j++)
		{
			HamiltonianMatrix(j,i) = BSplineObject->BSplineDerivative2OverlapIntegral(i,j);
			OverlapMatrix(j,i) = BSplineObject->BSplineOverlapIntegral(i,j);

			//Get lower triangle by transposing
			if (i != j)
			{
				HamiltonianMatrix(i,j) = HamiltonianMatrix(j,i);
				OverlapMatrix(i,j) = OverlapMatrix(j,i);
			}
		}
	}
	
	// Scale by mass (momentum squared term)
	HamiltonianMatrix *= -1.0 / (2.0 * Mass);
}


template<int Rank>
void Propagator<Rank>::ComputeHamiltonianEigenvectors(blitz::Array<cplx, 2> &HamiltonMatrix)
{
	using namespace blitz;
	
	// Resize Eigenvectors and Eigenvalues
	int N = BSplineObject->NumberOfBSplines;
	Eigenvalues.resize(N);
	Eigenvectors.resize(N,N);
	EigenvectorsInverse.resize(N,N);

	//Lapack stuff
	Array<int, 1> pivot(N);
	pivot = 0;
	linalg::LAPACK<cplx> lapack;

	//Compute LU factorization of overlap matrix
	lapack.CalculateLUFactorization(OverlapMatrix, pivot);

	//Invert overlap matrix
	lapack.CalculateMatrixInverse(OverlapMatrix, pivot);

	//Matrix-matrix multiply: inv(OverlapMat) * DiffMat 
	Array<cplx, 2> Matrix(N,N);
	Matrix = sum(OverlapMatrix(tensor::k, tensor::j) * HamiltonMatrix(tensor::i, tensor::k), tensor::k);

	// Store hamilton matrix. Matrix is in col-major order
	// (FORTRAN), so we transpose to get HamiltonianMatrix
	// in row-major order.
	HamiltonianMatrix.resize(N,N);
	HamiltonianMatrix = Matrix(tensor::j, tensor::i);

	//Calculate left and right eigenvectors of Matrix 
	//(inverse overlap times differentiation matrix) 
	bool calculateLeft = true;
	bool calculateRight = true;
	lapack.CalculateEigenvectorFactorization(calculateLeft, calculateRight, Matrix, Eigenvalues, EigenvectorsInverse, Eigenvectors);
}

//Instantiate propagators of rank 1-4
template class Propagator<1>;
template class Propagator<2>;
template class Propagator<3>;
template class Propagator<4>;

} //Namespace


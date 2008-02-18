#include "bsplinepropagator.h"
#include "../../representation/representation.h"
#include "../../representation/combinedrepresentation.h"
#include "../../representation/bspline/bsplinerepresentation.h"
#include "../../utility/blitzblas.h"
#include "../../utility/blitztricks.h"
#include "../../utility/blitzlapack.h"



namespace BSpline
{

cplx ToCplx(double a)
{
	return cplx(a, 0);
}

BZ_DECLARE_FUNCTION(ToCplx);


template<int Rank>
void Propagator<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	config.Get("mass", Mass);
	cout << "BSplinePropagator: Mass = " << Mass << endl;
}


template<int Rank>
void Propagator<Rank>::Setup(const cplx &dt, const Wavefunction<Rank> &psi, int rank)
{
	using namespace blitz;

	//Get b-spline object from wavefunction
	typename CombinedRepresentation<Rank>::Ptr repr = dynamic_pointer_cast< CombinedRepresentation<Rank> >(psi.GetRepresentation());
	BSplineRepresentation::Ptr subRepr = dynamic_pointer_cast< BSplineRepresentation >(repr->GetRepresentation(rank));
	BSplineObject = subRepr->GetBSplineObject();

	//Set class parameters
	PropagateRank = rank;

	int N = psi.GetRepresentation()->GetFullShape()(rank);

	//create some temporary arrays
	Array<cplx, 2> HamiltonianMatrixSetup;

	//Call setup routine to calculate Hamiltonian matrix
	SetupHamiltonianMatrix(HamiltonianMatrixSetup);

	//Obtain eigenvectors for HamiltonianMatrixSetup
	ComputeHamiltonianEigenvectors(HamiltonianMatrixSetup);

	//Create transform-and-multiply-matrix matrices
	PropagationMatrix.resize(Eigenvectors.shape());
	HamiltonianMatrix.resize(Eigenvectors.shape());

	cplx hamDot = 0;
	cplx expDot = 0;

	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			hamDot = 0;
			expDot = 0;
			for (int k=0; k<N; k++)
			{
				expDot += Eigenvectors(k,i) * exp(- I * dt * Eigenvalues(k)) * Eigenvectors(k,j);
				hamDot += Eigenvectors(k,i) * Eigenvalues(k) * Eigenvectors(k,j);
			}
			PropagationMatrix(i,j) = expDot;
			HamiltonianMatrix(i, j) = hamDot;
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

	//TODO: Make this faster, much faster, by calling
	//on some vectorized function in blas.
	
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
void Propagator<Rank>::SetupHamiltonianMatrix(blitz::Array<cplx, 2> HamiltonianMatrix)
{
	// I = ka+1+j
	// J = j
	
	int k = BSplineObject->MaxSplineOrder;
	int N = BSplineObject->NumberOfBSplines;
	int ldab = k + 1;
	HamiltonianMatrix.resize(N, ldab);

	for (int i = 0; i < N; i++)
	{
		int jMax = std::min(i+k, N);
		for (int j = i; j < jMax; j++)
		{
			int lapackI = k - 1 + i - j;
			int lapackJ = j;
			HamiltonianMatrix(lapackJ,lapackI) = BSplineObject->BSplineDerivative2OverlapIntegral(i,j);
		}
	}
	
	// Scale by mass (momentum squared term)
	HamiltonianMatrix *= -1.0 / (2.0 * Mass);
}

template<int Rank>
void Propagator<Rank>::ComputeHamiltonianEigenvectors(blitz::Array<cplx, 2> HamiltonMatrix)
{
	using namespace blitz;
	
	bool computeEigenvectors = true;
	linalg::LAPACK<cplx>::HermitianStorage upperOrLower = linalg::LAPACK<cplx>::HermitianUpper;
	int N = BSplineObject->NumberOfBSplines;

	// Resize Eigenvectors and Eigenvalues
	Eigenvalues.resize(N);
	Eigenvectors.resize(N,N);

	// Copy overlap matrix from BSplineObject, since LAPACK routine makes changes to it
	//Array<double, 2> doubleOverlapMatrix = BSplineObject->GetBSplineOverlapMatrix();
	//Array<cplx, 2> overlapMatrix(doubleOverlapMatrix.shape);
	//overlapMatrix = doubleOverlapMatrix(tensor::i);
	Array<cplx, 2> overlapMatrix ( BSplineObject->GetBSplineOverlapMatrix()(tensor::i) );
	
	// Solve generalized eigenvalue problem using LAPACK routine zhbgv
	linalg::LAPACK<cplx> lapack;
	lapack.CalculateEigenvectorFactorizationGeneralizedHermitianBanded(computeEigenvectors,
	                                                            upperOrLower,
	                                                            HamiltonMatrix,
	                                                            overlapMatrix,
	                                                            Eigenvalues,
                                                                Eigenvectors);
}


template class Propagator<1>;
template class Propagator<2>;
template class Propagator<3>;
template class Propagator<4>;


} //Namespace


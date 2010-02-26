#include "orthopolpropagator.h"
#include "../../utility/blitzblas.h"
#include "../../utility/blitztricks.h"
#include "../../representation/representation.h"

namespace OrthoPol
{
using namespace blitz;

cplx ToCplx(double a)
{
	return cplx(a, 0);
}

BZ_DECLARE_FUNCTION(ToCplx);


template<int Rank>
void Propagator<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	config.Get("mass", Mass);
	cout << "OrthoPolPropagator: Mass = " << Mass << endl;
}

template<int Rank>
void Propagator<Rank>::Setup(const Parameter &param, const cplx &dt, const Wavefunction<Rank> &psi, int rank)
{
	//Set class parameters
	N = psi.GetRepresentation()->GetFullShape()(rank);
	PropagateRank = rank;
	Param = param;
	
	// Call setup routines to create propagation matrix
	Array<double, 1> x;
	Array<double, 1> w;
	OrthoPolSetup(N, Param, Eigenvalues, Eigenvectors, x, w);

/*
	for (int i=0; i<N; i+=2)
	{
		Eigenvectors(i, Range::all()) = 0;
	}
*/

	PropagationMatrix.resize(N,N);
	PropagationMatrix = 0;
	DiffMatrix.resize(N,N);
	DiffMatrix = 0;
	cplx diffDot = 0;
	cplx expDot = 0;
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			diffDot = 0;
			expDot = 0;
			for (int k=0; k<N; k++)
			{
				expDot += Eigenvectors(k,i) * exp(- I * dt * Eigenvalues(k)) * Eigenvectors(k,j);
				diffDot += Eigenvectors(k,i) * Eigenvalues(k) * Eigenvectors(k,j);
			}
			PropagationMatrix(i,j) = expDot * w(j);
			DiffMatrix(i, j) = diffDot * w(j);
		}
	}

	// Allocate temp data
	TempData.resize(N);
}

template<int Rank>
void Propagator<Rank>::AdvanceStep(Wavefunction<Rank> &psi)
{
	//Map the data to a 3D array, where the radial part is the 
	//middle rank
	Array<cplx, 3> data3d = MapToRank3(psi.Data, PropagateRank, 1);

	//Propagate the 3D array
	ApplyPropagationMatrix(data3d);
}

template<int Rank>
void Propagator<Rank>::ApplyPropagationMatrix(Array<cplx, 3> &data)
{
	//TODO: Make this faster, much faster, by calling
	//on some vectorized function in blas.
	
	Array<cplx, 1> temp = TempData;
	Array<cplx, 2> prop = PropagationMatrix;
			
	//iterate over the array direction which is not propagated
	for (int i=0; i<data.extent(0); i++)
	{
		for (int j=0; j<data.extent(2); j++)
		{
			Array<cplx, 1> v = data(i, Range::all(), j);
			temp = v; // v.copy();
			MatrixVectorMultiply(prop, temp, v);
		}
	}
}

template<int Rank>
void Propagator<Rank>::ApplyDifferentiationMatrix(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi)
{
	//Map the data to a 3D array, where the radial part is the 
	//middle rank
	Array<cplx, 3> srcData = MapToRank3(srcPsi.Data, PropagateRank, 1);
	Array<cplx, 3> dstData = MapToRank3(dstPsi.Data, PropagateRank, 1);
	Array<cplx, 2> matrix = GetDifferentiationMatrix();

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

template class Propagator<1>;
template class Propagator<2>;
template class Propagator<3>;
template class Propagator<4>;

};


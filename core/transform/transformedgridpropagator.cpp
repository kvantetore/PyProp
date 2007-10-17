#include "transformedgridpropagator.h"
#include "../utility/blitzblas.h"
#include "../utility/blitztricks.h"
#include "../representation/representation.h"

namespace TransformedGrid
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
	cout << "TransformedGridPropagator: Mass = " << Mass << endl;
}

template<int Rank>
void Propagator<Rank>::Setup(const Parameter &param, const cplx &dt, const Wavefunction<Rank> &psi, int rank)
{
	firstIndex i;
	secondIndex j;
	thirdIndex k;

	//Set class parameters
	N = psi.GetRepresentation()->GetFullShape()(rank);
	PropagateRank = rank;
	Param = param;

	//create some temporary arrays
	Array<double, 1> ewReal;    //eigenvalue 
	Array<double, 2> evReal;    //eigenvector matrix
	Array<double, 2> evInvReal; //inverse eigenvector matrix

	//Call setup routines to get the eigenvector decomposition
	setup(N, param, ewReal, evReal, evInvReal);

	//We need the eigenvalues and vectors in complex format
	Array<cplx, 1> ew(ewReal.shape());    //eigenvalue 
	Array<cplx, 2> ev(evReal.shape());    //eigenvector matrix
	Array<cplx, 2> evExp(evReal.shape());     //eigenvectors scaled by exp(ew)
	Array<cplx, 2> evDiff(evReal.shape());    //eigenvectors scaled by ew
	Array<cplx, 2> evInv(evInvReal.shape());  //inverse eigenvector matrix

	ew = ewReal(tensor::i);
	ev = evReal(tensor::i, tensor::j);
	evInv = evInvReal(tensor::i, tensor::j);

	//scale eigenvectors by complex rotation
	//The missing minus sign in the exponent is included in the matrix.
	evExp = ev(i,j) * exp( I * dt * ew(j) / (2.0 * Mass));
	evDiff = - ev(i,j) * ew(j) / (2.0 * Mass);

	//Create full matrix to propagate wavefunction
	PropagationMatrix.resize(ev.shape());
	MatrixMatrixMultiply(evExp, evInv, PropagationMatrix);

	//Create full differentiation matrix
	DiffMatrix.resize(ev.shape());
	MatrixMatrixMultiply(evDiff, evInv, DiffMatrix);

	//Allocate temp data
	TempData.resize(ew.extent(0));
}

template<int Rank>
void Propagator<Rank>::AdvanceStep(Wavefunction<Rank> &psi)
{
	//Map the data to a 3D array, where the radial part is the 
	//middle rank
	Array<cplx, 3> data3d = MapToRank3(psi.Data, PropagateRank, 1);

	//Propagate the 3D array
	ApplyMatrix(PropagationMatrix, data3d);
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

template<int Rank>
void Propagator<Rank>::ApplyMatrix(const Array<cplx, 2> &matrix, Array<cplx, 3> &data)
{
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

template<int Rank>
blitz::Array<cplx, 2> Propagator<Rank>::GetDifferentiationMatrix()
{
	//TODO: Move initialization of DiffMatrix in here to save memory in case
	//we don't need it.
	return DiffMatrix;
}



template class Propagator<1>;
template class Propagator<2>;
template class Propagator<3>;
template class Propagator<4>;

};


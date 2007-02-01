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
void Propagator<Rank>::Setup(const Parameter &param, const cplx &dt, const Wavefunction<Rank> &psi, int rank)
{
	firstIndex i;
	secondIndex j;
	thirdIndex k;

	//Set class parameters
	N = psi.GetRepresentation().GetFullShape()(rank);
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
	Array<cplx, 2> evInv(evInvReal.shape());  //inverse eigenvector matrix

			
	ew = ewReal(tensor::i);
	ev = evReal(tensor::i, tensor::j);
	evInv = evInvReal(tensor::i, tensor::j);

//	cout.precision(15);
//	cout << "Eigenvalues: " << ewReal << endl;
//	cout << "Eigenvectors: " << evReal << endl;
//	cout << "Eigenvectorsinv: " << evInvReal << endl;
	
	
	//scale eigenvectors by complex rotation
	//The missing minus sign in the exponent is included in the matrix.
	double mass = 1.0;
	ev = ev(i,j) * exp( I * dt * ew(j) / (2.0 * mass));

	//Create full matrix to propagate wavefunction
	PropagationMatrix.resize(ev.shape());
	MatrixMatrixMultiply(ev, evInv, PropagationMatrix);

	//Allocate temp data
	TempData.resize(ew.extent(0));
}

template<int Rank>
void Propagator<Rank>::AdvanceStep(Wavefunction<Rank> &psi)
{
	if (PropagateRank != 0)
	{
		cout << "Currently TransformedGridPropagator only works for the first rank" <<endl;
		throw std::runtime_error("Invalid rank");
	}

	//Map the data to a 2D array, where the radial part is the 
	//first rank (of highest stride)
	Array<cplx, 2> data2d = MapToRank2(psi.Data, 1);

	//Propagate the 2D array
	ApplyPropagationMatrix(data2d);
}

template<int Rank>
void Propagator<Rank>::ApplyPropagationMatrix(Array<cplx, 2> data)
{
	//TODO: Make this faster, much faster, by calling
	//on some vectorized function in blas.
	
	Array<cplx, 1> temp = TempData;
	Array<cplx, 2> prop = PropagationMatrix;
			
	//iterate over the array direction which is not propagated
	for (int i=0; i<data.extent(1); i++)
	{
		Array<cplx, 1> v = data(Range::all(), i);
		temp = v.copy();
		MatrixVectorMultiply(prop, temp, v);
	}
}


template class Propagator<1>;
template class Propagator<2>;
template class Propagator<3>;
template class Propagator<4>;
template class Propagator<5>;
template class Propagator<6>;


};


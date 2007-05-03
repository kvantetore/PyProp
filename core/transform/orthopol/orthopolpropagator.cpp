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
	N = psi.GetRepresentation().GetFullShape()(rank);
	PropagateRank = rank;
	Param = param;
	
	// Call setup routines to create propagation matrix
	OrthoPolSetup(N, Param, dt, PropagationMatrix);

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


template class Propagator<1>;
template class Propagator<2>;
template class Propagator<3>;
template class Propagator<4>;

};


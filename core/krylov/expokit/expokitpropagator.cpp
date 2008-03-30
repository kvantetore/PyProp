#include "expokitpropagator.h"
#include "f2c/expokit.h"

#include <core/wavefunction.h>

#include "../krylovcommon.h"

namespace krylov
{

using namespace boost::python;


template<int Rank>
void ExpokitPropagator<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	config.Get("krylov_basis_size", BasisSize);
	config.Get("krylov_tolerance", Tolerance);
	config.Get("krylov_norm", MatrixNorm);
}


template<int Rank>
void ExpokitPropagator<Rank>::Setup(const Wavefunction<Rank> &psi)
{
	int n = psi.GetData().size();
	int m = BasisSize;
	int workspaceSize = n * (m+2) + 5 *(m+2)*(m+2) + 6 + 1;
	cout << "Allocating krylov workspace of " << (double) workspaceSize * sizeof(cplx) / (1024.0*1024.0) << "MB" << endl;
	this->Workspace.resize(workspaceSize);
	this->IntegerWorkspace.resize(m+2);
}


template<int Rank>
void ExpokitPropagator<Rank>::AdvanceStep(object callback, Wavefunction<Rank> &psi, Wavefunction<Rank> &tempPsi, cplx dt, double t)
{
	int trace = 0; //Set to 1 to print diagnostics while propagating
	int flag = 0;  //Should be 0 on return

	int n = psi.GetData().size();	
	int m = BasisSize;

	//in, out
	cplx* in = psi.GetData().data();
	cplx* out = tempPsi.GetData().data();

	double inNorm = psi.GetNorm();
	//cout << "inNorm = " << inNorm << endl;
	if (inNorm < 1e-6)
	{
		cout << "Wavefunction is zero: " << endl;
		cout << "Initial norm = " << inNorm << endl;
	}

	//Set up class variables needed for callback
	this->Psi = &psi;
	this->TempPsi = &tempPsi;
	this->MultiplyCallback = callback;
	//Set up timestep
	this->TimeStep = sqrt(sqr(imag(dt) + sqr(real(dt))));
	//cout << "Using timestep " << this->TimeStep << endl;
	this->ImaginaryTime = false;
	if (abs(imag(dt)) > 1e-10)
	{
		this->ImaginaryTime = true;
	}
	this->CurTime = t;

	int wspSize = Workspace.size();
	int iwspSize = IntegerWorkspace.size();

	expokit::zgexpv( 
		&n, 
		&m, 
		&this->TimeStep, 
		in, 
		out, 
		&Tolerance, 
		&MatrixNorm, 
		Workspace.data(), 
		&wspSize, 
		IntegerWorkspace.data(), 
		&iwspSize, 
		MultiplyHamiltonian<Rank>, 
		this, 
		&trace, 
		&flag
	);

	if (flag != 0)
	{
		cout << "Error: Expokit retured nonzero value " << flag << endl;
		throw std::runtime_error("Expokit returned nonzero value");
	}

	//Copy output data back to wavefunction
	psi.GetData() = tempPsi.GetData();
}

template class ExpokitPropagator<1>;
template class ExpokitPropagator<2>;
template class ExpokitPropagator<3>;
template class ExpokitPropagator<4>;

} //Namespace
	

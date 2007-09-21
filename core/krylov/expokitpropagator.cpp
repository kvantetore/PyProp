#include "expokitpropagator.h"
#include "expokit/expokit.h"

#include <core/wavefunction.h>

namespace krylov
{

using namespace boost::python;


template<int Rank>
void MultiplyHamiltonian(void *data, cplx *inBuffer, cplx *outBuffer)
{
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;

	if (inBuffer == outBuffer)
	{
		cout << "Error: both input and output-buffer is the same." << endl;
		cout << "  iBuffer address = " << inBuffer << endl;
		cout << "  oBuffer address = " << outBuffer << endl;
	}

	//Some local variables for simplicity
	ExpokitPropagator<Rank> *propagator = static_cast<ExpokitPropagator<Rank>*>(data);
	Wavefunction<Rank> *psi = propagator->Psi;
	Wavefunction<Rank> *tempPsi = propagator->TempPsi;

	//Wrap the data buffers in blitz arrays and initialize outdata to 0
	DataVector shape = psi->GetData().shape();
	DataVector stride = psi->GetData().stride();
	DataArray inData(inBuffer, shape, stride, blitz::neverDeleteData);
	DataArray outData(outBuffer, shape, stride, blitz::neverDeleteData);
	outData = 0;

	cplx inNorm = sum((inData * conj(inData)));
	//cout << "In Norm = " << inNorm << endl;

	//Remember the original data buffers of psi and tempPsi
	DataArray psiOrigData(psi->GetData());
	DataArray tempPsiOrigData(tempPsi->GetData());
	//Use the data buffers supplied by expokit
	psi->SetData(inData);
	tempPsi->SetData(outData);

	//Perform Hamilton-Wavefunction multiplication
	propagator->MultiplyCallback(psi, tempPsi, propagator->TimeStep, propagator->CurTime);

	//Restore the databuffers of psi and tempPsi
	psi->SetData(psiOrigData);
	tempPsi->SetData(tempPsiOrigData);

	//Scale with negative imaginary unit if we do NOT use imaginary time.
	if (!propagator->ImaginaryTime)
	{
		outData *= cplx(0.0, -1.0);
	}
	else
	{
		outData *= -1.0;
	}

	cplx outNorm = sum((outData * conj(outData)));
	if (real(outNorm) < 0.1)
	{
		cout << "Outnorm is very small. Something is amiss." << endl;
		cout << "  InNorm = " << inNorm << endl;
		cout << " outNorm = " << outNorm << endl;
		cout << "  inAddr = " << inBuffer << endl;
		cout << " outAddr = " << outBuffer << endl;
		cout << "    diff = " << (int)(inBuffer - outBuffer) << endl;
	}
}


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
	cout << "inNorm = " << inNorm << endl;
	if (inNorm < 1e-6)
	{
		cout << "Wavefunction is zero: " << endl;
		cout << "Initial norm = " << inNorm << endl;
	}

	//Set up class variables needed for callback
	Psi = &psi;
	TempPsi = &tempPsi;
	MultiplyCallback = callback;
	//Set up timestep
	TimeStep = sqrt(sqr(imag(dt) + sqr(real(dt))));
	cout << "Using timestep " << TimeStep << endl;
	ImaginaryTime = false;
	if (abs(imag(dt)) > 1e-10)
	{
		ImaginaryTime = true;
	}
	CurTime = t;

	int wspSize = Workspace.size();
	int iwspSize = IntegerWorkspace.size();

	expokit::zgexpv( 
		&n, 
		&m, 
		&TimeStep, 
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
	

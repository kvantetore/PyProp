#include "odewrapper.h"

#include <core/wavefunction.h>
#include <core/utility/fortran.h>

namespace ODE
{

typedef void (*callbackFunction)(double *t, cplx* in, cplx* out, void* data);

extern "C"
{
	void FORTRAN_NAME(ode)(callbackFunction callback, void* data, int* n, cplx* in, double* t, double* tout, double* relerr, double* abserr, int* flag, cplx* work, int* iwork);
}


// "Thin" wrapper of ode integrator
void cOde(callbackFunction callback, void* data, int n, cplx* in, double &t, double tout, double relerr, double abserr, int &flag, cplx* work, int* iwork )
{
	FORTRAN_NAME(ode)(callback, data, &n, in, &t, &tout, &relerr, &abserr, &flag, work, iwork);
}

using namespace boost::python;

template<int Rank>
void MultiplyHamiltonian(double *t, cplx *inBuffer, cplx *outBuffer, void *data)
{
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;

	//Some local variables for simplicity
	OdeWrapper<Rank> *propagator = static_cast< OdeWrapper<Rank>* >(data);
	Wavefunction<Rank> *psi = propagator->Psi;
	Wavefunction<Rank> *tempPsi = propagator->TempPsi;

	//Wrap the data buffers in blitz arrays and initialize outdata to 0
	DataVector shape = psi->GetData().shape();;
	DataVector stride = psi->GetData().stride();
	DataArray inData(inBuffer, shape, stride, blitz::neverDeleteData);
	DataArray outData(outBuffer, shape, stride, blitz::neverDeleteData);
	outData = 0;

	//Remember the original data buffers of psi and tempPsi
	DataArray psiOrigData(psi->GetData());
	DataArray tempPsiOrigData(tempPsi->GetData());

	//Use the data buffers supplied by expokit
	psi->SetData(inData);
	tempPsi->SetData(outData);

	//Perform Hamilton-Wavefunction multiplication
	object callback = propagator->MultiplyCallback;
	callback(psi, tempPsi, *t);

	//Restore the databuffers of psi and tempPsi
	psi->SetData(psiOrigData);
	tempPsi->SetData(tempPsiOrigData);

	if (propagator->ImTime)
	{
		outData = -outData;
	}
	else
	{
		//Scale with negative imaginary unit since this is not done by matrix-vector routine
		outData *= cplx(0.0, -1.0);
	}

	//cout << "Output2 = " << outData << endl;
	//throw std::runtime_error("STOP!");
}


template<int Rank>
void OdeWrapper<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	this->RelativeError = 1e-9;
	if (config.HasValue("relative_error"))
	{
		config.Get("relative_error", this->RelativeError);
		cout << "Using Relative error " << this->RelativeError << endl;
	}
	else
	{
		cout << "Using Default Relative error " << this->RelativeError << endl;
	}

	this->AbsoluteError = 1e-9;
	if (config.HasValue("absolute_error"))
	{
		config.Get("absolute_error", this->AbsoluteError);
		cout << "Using Absolute error " << this->AbsoluteError << endl;
	}
	else
	{
		cout << "Using Default Absolute error " << this->AbsoluteError << endl;
	}
}


template<int Rank>
void OdeWrapper<Rank>::Setup(const Wavefunction<Rank> &psi)
{
	int n = psi.GetData().size();
	int workspaceSize = 100 + 21 * n;
	cout << "Allocating ODE workspace of " << (double) workspaceSize * sizeof(cplx) / (1024.0*1024.0) << " MB" << endl;
	this->Work.resize(workspaceSize);
	this->Iwork.resize(50);
	this->Flag = 1;
	this->OutputTime = 0;
}


template<int Rank>
void OdeWrapper<Rank>::AdvanceStep(object callback, Wavefunction<Rank> &psi, Wavefunction<Rank> &tempPsi, cplx dt, double t)
{
	int n = 2 * psi.GetData().size();	

	//in, out
	cplx* in = psi.GetData().data();

	//Set up class variables needed for callback
	this->Psi = &psi;
	this->TempPsi = &tempPsi;
	this->MultiplyCallback = callback;

	//double tout = t + sqrt(sqr(imag(dt) + sqr(real(dt))));
	this->ImTime = std::abs(imag(dt)) > 1e-10;
	double tout = t + std::abs(dt);

	cOde(MultiplyHamiltonian<Rank>, this, n, in, this->OutputTime, tout, this->RelativeError, this->AbsoluteError, this->Flag, this->Work.data(), this->Iwork.data());

	if (this->Flag != 2)
	{
		cout << "Fatal error!" << endl;
		cout << "ODE returned FLAG = " << this->Flag << endl;
		throw std::runtime_error("ODE returned an error");
	}

	//Remove reference to callback to prevent locking memory on exit
	this->MultiplyCallback = object();
}

template class OdeWrapper<1>;
template class OdeWrapper<2>;
} //Namespace


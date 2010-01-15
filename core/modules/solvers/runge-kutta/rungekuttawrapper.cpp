#include "rungekuttawrapper.h"
#include <gsl/gsl_errno.h>
#include <core/wavefunction.h>

namespace RungeKutta
{

using namespace boost::python;

template<int Rank>
int MultiplyHamiltonian(double t, const double *inBuffer, double *outBuffer, void *data)
{
	typedef blitz::TinyVector<int, Rank> DataVector;
	typedef blitz::Array<cplx, Rank> DataArray;

	//Some local variables for simplicity
	RungeKuttaWrapper<Rank> *propagator = static_cast< RungeKuttaWrapper<Rank>* >(data);
	typename Wavefunction<Rank>::Ptr psi = propagator->Psi;
	typename Wavefunction<Rank>::Ptr tempPsi = propagator->TempPsi;

	//Wrap the data buffers in blitz arrays and initialize outdata to 0
	DataVector shape = psi->GetData().shape();;
	DataVector stride = psi->GetData().stride();
	DataArray inData((cplx *)inBuffer, shape, stride, blitz::neverDeleteData);
	DataArray outData((cplx *)outBuffer, shape, stride, blitz::neverDeleteData);
	outData = 0;

	//Remember the original data buffers of psi and tempPsi
	DataArray psiOrigData(psi->GetData());
	DataArray tempPsiOrigData(tempPsi->GetData());

	//Use the data buffers supplied by expokit
	psi->SetData(inData);
	tempPsi->SetData(outData);

	//Perform Hamilton-Wavefunction multiplication
	object callback = propagator->MultiplyCallback;
	callback(psi, tempPsi, t);

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

	return GSL_SUCCESS;
}


template<int Rank>
void RungeKuttaWrapper<Rank>::ApplyConfigSection(const ConfigSection &config)
{

	config.Get("integrator_type", this->IntegratorType);

}


template<int Rank>
void RungeKuttaWrapper<Rank>::Setup(typename Wavefunction<Rank>::Ptr psi)
{

	size_t n = psi->GetData().size() * 2;

	//Set integrator type
	if (IntegratorType == IntegratorRK2) { IntegratorTypeObject = gsl_odeiv_step_rk2; }
	else if (IntegratorType == IntegratorRK4) { IntegratorTypeObject = gsl_odeiv_step_rk4; }
	else if (IntegratorType == IntegratorRKF45) { IntegratorTypeObject = gsl_odeiv_step_rkf45; }
	else if (IntegratorType == IntegratorRKCK) { IntegratorTypeObject = gsl_odeiv_step_rkck; }
	else if (IntegratorType == IntegratorRK8PD) { IntegratorTypeObject = gsl_odeiv_step_rk8pd; }

	else
	{
		cout << "Unknown integration algorithm: " << IntegratorType;
	}

	//Allocate integrator workspace
	Integrator = gsl_odeiv_step_alloc (IntegratorTypeObject, n);
	cout << "Using integrator: " << gsl_odeiv_step_name(Integrator) << endl;

	//Setup integrator data struct
	sys.function = MultiplyHamiltonian<Rank>;
	sys.jacobian = 0;
	sys.dimension = n;
	sys.params = this;
}


template<int Rank>
void RungeKuttaWrapper<Rank>::AdvanceStep(object callback, typename Wavefunction<Rank>::Ptr psi, typename Wavefunction<Rank>::Ptr tempPsi, cplx dt, double t)
{
	//in/out and element absolute error
	double* inout = (double *)psi->GetData().data();
	double* y_err = (double *)tempPsi->GetData().data();

	//Set up class variables needed for callback
	this->Psi = psi;
	this->TempPsi = tempPsi;
	this->MultiplyCallback = callback;

	double timeStep;
	if (std::abs(imag(dt)) > 1e-10)
	{
		this->ImTime = true;
		timeStep = std::abs(dt);
	}
	else
	{
		this->ImTime = false;
		timeStep = std::real(dt);
	}

	int status = gsl_odeiv_step_apply(Integrator, t, timeStep, inout, y_err, 0, 0, &sys);

	if (status != GSL_SUCCESS)
	{
		cout << "Integrator returned an error: " << gsl_strerror (status) << endl;
	}

	//Remove reference to callback to prevent locking memory on exit
	this->MultiplyCallback = object();
}

template class RungeKuttaWrapper<1>;
template class RungeKuttaWrapper<2>;
} //Namespace


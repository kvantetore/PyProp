#include "pampwrapper.h"

#include <core/wavefunction.h>
#include <core/representation/representation.h>
#include <core/mpi/distributedmodel.h>
#include <core/utility/fortran.h>
#include <core/utility/blitztricks.h>

#include "../pypropfunctor.h"

namespace krylov
{

using namespace boost::python;

/* Implementation of PampWrapper */
template<int Rank>
void PampWrapper<Rank>::ApplyOperator(blitz::Array<cplx, 1> &input, blitz::Array<cplx, 1> &output)
{
	if (Psi == 0)
	{
		throw std::runtime_error("Psi is 0");
	}
	if (TempPsi == 0)
	{
		throw std::runtime_error("TempPsi is 0");
	}

	//Map the 1d vectors to a a blitz array of correct shape
	DataVector shape = Psi->GetData().shape();
	DataVector stride = Psi->GetData().stride();
	DataArray inData(input.data(), shape, stride, blitz::neverDeleteData);
	DataArray outData(output.data(), shape, stride, blitz::neverDeleteData);
	outData = 0;

	//Set psi and tempPsi to point to correct vectorsbuffers
	DataArray oldData = Psi->GetData();
	DataArray oldTempData = TempPsi->GetData();
	Psi->SetData(inData);
	TempPsi->SetData(outData);

	Callback(Psi, TempPsi, TimeStep, CurTime);

	//Restore the former buffers
	Psi->SetData(oldData);
	TempPsi->SetData(oldTempData);
}


template<int Rank>
void PampWrapper<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	config.Get("krylov_basis_size", Propagator.BasisSize);

	if (config.HasValue("krylov_exponentiation_method"))
	{
		int exponentiation;
		config.Get("krylov_exponentiation_method", exponentiation);
		Propagator.Exponentiation = (pamp::pAMP<cplx>::ExponentiationMethod) exponentiation;
	}
}


template<int Rank>
void PampWrapper<Rank>::Setup(const typename Wavefunction<Rank>::Ptr psi)
{
	Propagator.MatrixSize = psi->GetData().size();
    Propagator.DisableMPI = psi->GetRepresentation()->GetDistributedModel()->IsSingleProc();
	Propagator.MatrixOperator = typename PypropOperatorFunctor<Rank>::Ptr( new PypropOperatorFunctor<Rank>(this) );
	Propagator.Setup();
}


template<int Rank>
void PampWrapper<Rank>::AdvanceStep(object callback, typename Wavefunction<Rank>::Ptr psi, typename Wavefunction<Rank>::Ptr tempPsi, cplx dt, double t, bool usePypropIntegration)
{
	//The callback-functions uses these variables
	this->Psi = psi;
	this->TempPsi = tempPsi;
	this->Callback = callback;
	this->CurTime = t;
	this->TimeStep = dt;

	//Use our custom integration
	if (usePypropIntegration)
	{
		typename piram::IntegrationFunctor<cplx, double>::Ptr integration(new PypropIntegrationFunctor<Rank>(Psi, TempPsi));
		Propagator.Integration = integration;
	}

	typename Wavefunction<Rank>::DataArray vector = psi->GetData();
	blitz::Array<cplx, 1> vector1d = MapToRank1(vector);
	Propagator.PropagateVector(vector1d, dt);

	//Zero the pointers to avoid mishaps
	this->Psi = typename Wavefunction<Rank>::Ptr();
	this->TempPsi = typename Wavefunction<Rank>::Ptr();
}

template class PampWrapper<1>;
template class PampWrapper<2>;
template class PampWrapper<3>;
template class PampWrapper<4>;

} //Namespace
	

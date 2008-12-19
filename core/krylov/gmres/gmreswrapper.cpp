#include "gmreswrapper.h"

#include <core/wavefunction.h>
#include <core/representation/representation.h>
#include <core/mpi/distributedmodel.h>
#include <core/utility/fortran.h>
#include <core/utility/blitztricks.h>

#include "../pypropfunctor.h"

namespace krylov
{

using namespace boost::python;

/* Implementation of GmresWrapper */
template<int Rank>
void GmresWrapper<Rank>::ApplyOperator(blitz::Array<cplx, 1> &input, blitz::Array<cplx, 1> &output)
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
	DataArray inData2(input.data(), shape, stride, blitz::neverDeleteData);
	DataArray outData(output.data(), shape, stride, blitz::neverDeleteData);
	DataArray outData2(output.data(), shape, stride, blitz::neverDeleteData);
	outData = 0;

	//Set psi and tempPsi to point to correct vectorsbuffers
	DataArray oldData = Psi->GetData();
	DataArray oldTempData = TempPsi->GetData();
	Psi->SetData(inData);
	TempPsi->SetData(outData);

	OperatorCallback(Psi, TempPsi);

	//Restore the former buffers
	Psi->SetData(oldData);
	TempPsi->SetData(oldTempData);
}


template<int Rank>
void GmresWrapper<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	config.Get("krylov_basis_size", Solver.BasisSize);

	if (config.HasValue("krylov_tolerance"))
	{
		config.Get("krylov_tolerance", Solver.Tolerance);
	}

	//Perform double orthogonalization step?
	if (config.HasValue("krylov_double_orthogonalization"))
	{
		bool performDoubleOrthogonalization;
		config.Get("krylov_double_orthogonalization", performDoubleOrthogonalization);
		cout << "Using doubleorth = " << performDoubleOrthogonalization << endl;
		Solver.PerformDoubleOrthogonalization = performDoubleOrthogonalization;
	}
}


template<int Rank>
void GmresWrapper<Rank>::Setup(const typename Wavefunction<Rank>::Ptr psi)
{
	Solver.MatrixSize = psi->GetData().size();
    Solver.DisableMPI = psi->GetRepresentation()->GetDistributedModel()->IsSingleProc();
	Solver.MatrixOperator = typename PypropOperatorFunctor<Rank>::Ptr( new PypropOperatorFunctor<Rank>(this) );
	Solver.Setup();
}


template<int Rank>
void GmresWrapper<Rank>::Solve(object callback, typename Wavefunction<Rank>::Ptr psi, typename Wavefunction<Rank>::Ptr tempPsi, bool usePypropIntegration)
{
	//The callback-functions uses these variables
	this->Psi = psi;
	this->TempPsi = tempPsi;
	this->OperatorCallback = callback;

	//Use the inner product from pyprop, or assume that
	//the matrix is scaled appropriately
	if (usePypropIntegration)
	{
		typename piram::IntegrationFunctor<cplx, double>::Ptr integration(new PypropIntegrationFunctor<Rank>(Psi, TempPsi));
		Solver.Integration = integration;
	}

	typename Wavefunction<Rank>::DataArray inputVector = psi->GetData();
	blitz::Array<cplx, 1> inputVector1d = MapToRank1(inputVector);
	typename Wavefunction<Rank>::DataArray outputVector = tempPsi->GetData();
	blitz::Array<cplx, 1> outputVector1d = MapToRank1(outputVector);

	Solver.SolveVector(inputVector1d, outputVector1d);

	//Zero the pointers to avoid mishaps
	this->Psi = typename Wavefunction<Rank>::Ptr();
	this->TempPsi = typename Wavefunction<Rank>::Ptr();
}

template class GmresWrapper<1>;
template class GmresWrapper<2>;
template class GmresWrapper<3>;
template class GmresWrapper<4>;

} //Namespace
	

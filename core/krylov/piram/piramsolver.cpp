#include "piramsolver.h"

#include <core/wavefunction.h>
#include <core/representation/representation.h>
#include <core/mpi/distributedmodel.h>
#include <core/utility/fortran.h>
#include <core/utility/blitztricks.h>

namespace krylov
{

/* Implementation of PiramSolver */

template<int Rank>
void PiramSolver<Rank>::SetupResidual(blitz::Array<cplx, 1> &residual)
{
	if (Psi == 0)
	{
		throw std::runtime_error("Psi is 0");
	}

	residual = MapToRank1(Psi->GetData());
}

template<int Rank>
void PiramSolver<Rank>::ApplyOperator(blitz::Array<cplx, 1> &input, blitz::Array<cplx, 1> &output)
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

	Callback(Psi, TempPsi);

	//Restore the former buffers
	Psi->SetData(oldData);
	TempPsi->SetData(oldTempData);
}


template<int Rank>
void PiramSolver<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	config.Get("krylov_basis_size", Solver.BasisSize);
	config.Get("krylov_tolerance", Solver.Tolerance);
	config.Get("krylov_eigenvalue_count", Solver.EigenvalueCount);
	config.Get("krylov_max_iteration_count", Solver.MaxRestartCount);
	config.Get("krylov_use_random_start", Solver.UseRandomStart);

	if (config.HasValue("krylov_eigenvalue_shift"))
	{
		cplx shift;
		config.Get("krylov_eigenvalue_shift", shift);
	
		typedef piram::CompareComplexNearShift compareType;
		typedef piram::ShiftFunctorSelectRitzValues< cplx, compareType > shiftFunctorType;
		compareType compare(shift);
		Solver.CalculateShifts = typename shiftFunctorType::Ptr( new shiftFunctorType(compare) );
	}

	if (config.HasValue("inverse_iterations"))
	{
		cout << "PIRAM: Using inverse iterations" << endl;
		bool inverseIt;
		config.Get("inverse_iterations", inverseIt);

		if (inverseIt)
		{
			typedef piram::CompareComplexGreaterAbs compareType;
			typedef piram::ShiftFunctorSelectRitzValues< cplx, compareType > shiftFunctorType;
			compareType compare;
			Solver.CalculateShifts = typename shiftFunctorType::Ptr( new shiftFunctorType(compare) );
		}
	}
}


template<int Rank>
void PiramSolver<Rank>::Setup(const typename Wavefunction<Rank>::Ptr psi)
{
	Solver.MatrixSize = psi->GetData().size();
    Solver.DisableMPI = psi->GetRepresentation()->GetDistributedModel()->IsSingleProc();

	Solver.SetupResidual = typename PypropSetupResidualFunctor<Rank>::Ptr( new PypropSetupResidualFunctor<Rank>(this) );
	Solver.MatrixOperator = typename PypropOperatorFunctor<Rank>::Ptr( new PypropOperatorFunctor<Rank>(this) );

	Solver.Setup();
}


template<int Rank>
void PiramSolver<Rank>::Solve(object callback, typename Wavefunction<Rank>::Ptr psi, typename Wavefunction<Rank>::Ptr tempPsi, bool usePypropIntegration)
{
	//The callback-functions uses these variables
	this->Psi = psi;
	this->TempPsi = tempPsi;
	this->Callback = callback;

	//Use our custom integration
	if (usePypropIntegration)
	{
		typename piram::IntegrationFunctor<cplx, double>::Ptr integration(new PypropIntegrationFunctor<Rank>(Psi, TempPsi));
		Solver.Integration = integration;
	}
	
	Solver.Solve();
	Solver.Postprocess();

	//Zero the pointers to avoid mishaps
	this->Psi = typename Wavefunction<Rank>::Ptr();
	this->TempPsi = typename Wavefunction<Rank>::Ptr();
}

template class PiramSolver<1>;
template class PiramSolver<2>;
template class PiramSolver<3>;
template class PiramSolver<4>;

} //Namespace
	

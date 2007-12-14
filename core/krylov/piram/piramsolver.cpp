#include "piramsolver.h"

#include <core/wavefunction.h>
#include <core/representation/representation.h>
#include <core/mpi/distributedmodel.h>
#include <core/utility/fortran.h>
#include <core/utility/blitztricks.h>

namespace krylov
{

using namespace boost::python;

/* Functors for pIRAM */
template<int Rank>
class PypropOperatorFunctor : public piram::OperatorFunctor<cplx>
{
public:
	typedef boost::shared_ptr< PypropOperatorFunctor > Ptr;

	PypropOperatorFunctor(PiramSolver<Rank> *solver) : Solver(solver) {}

	virtual void operator()(blitz::Array<cplx, 1> &in, blitz::Array<cplx, 1> &out) 
	{
		Solver->ApplyOperator(in, out);
	}

private:
	PiramSolver<Rank>* Solver;
};

template<int Rank>
class PypropSetupResidualFunctor : public piram::SetupResidualFunctor<cplx>
{
public:
	typedef boost::shared_ptr< PypropSetupResidualFunctor > Ptr;

	PypropSetupResidualFunctor(PiramSolver<Rank> *solver) : Solver(solver) {}

	virtual void operator()(blitz::Array<cplx, 1> &residual) 
	{
		Solver->SetupResidual(residual);
	}

private:
	PiramSolver<Rank>* Solver;
};


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
	Psi->SetData(inData);
	TempPsi->SetData(outData);

	Callback(Psi, TempPsi);
}


template<int Rank>
void PiramSolver<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	config.Get("krylov_basis_size", Solver.BasisSize);
	config.Get("krylov_tolerance", Solver.Tolerance);
	config.Get("krylov_eigenvalue_count", Solver.EigenvalueCount);
	config.Get("krylov_max_iteration_count", Solver.MaxRestartCount);
	config.Get("krylov_use_random_start", Solver.UseRandomStart);
}


template<int Rank>
void PiramSolver<Rank>::Setup(const Wavefunction<Rank> &psi)
{
	Solver.MatrixSize = psi.GetData().size();
    Solver.DisableMPI = psi.GetRepresentation()->GetDistributedModel()->IsSingleProc();

	Solver.SetupResidual = typename PypropSetupResidualFunctor<Rank>::Ptr( new PypropSetupResidualFunctor<Rank>(this) );
	Solver.MatrixOperator = typename PypropOperatorFunctor<Rank>::Ptr( new PypropOperatorFunctor<Rank>(this) );

	Solver.Setup();
}


template<int Rank>
void PiramSolver<Rank>::Solve(object callback, Wavefunction<Rank> &psi, Wavefunction<Rank> &tempPsi)
{
	//The callback-functions uses these variables
	this->Psi = &psi;
	this->TempPsi = &tempPsi;
	this->Callback = callback;

	Solver.Solve();
	Solver.Postprocess();

	//Zero the pointers to avoid mishaps
	this->Psi = 0;
	this->TempPsi = 0;
}

template class PiramSolver<1>;
template class PiramSolver<2>;
template class PiramSolver<3>;
template class PiramSolver<4>;

} //Namespace
	

#ifndef BSPLINEPROPAGATOR_H
#define BSPLINEPROPAGATOR_H

#include "../../common.h"
#include "../../wavefunction.h"
#include "bspline.h"

namespace BSpline
{

template<int Rank>
class Propagator
{
private:
	BSpline::Ptr BSplineObject;

	blitz::Array<cplx, 2> PropagationMatrix;     // LAPACK matrix
	blitz::Array<cplx, 2> PropagationMatrixBlas; // BLAS matrix
	blitz::Array<cplx, 2> HamiltonianMatrix;     // BLAS matrix
	blitz::Array<cplx, 2> OverlapMatrix;         // LAPACK matrix
	blitz::Array<cplx, 1> TempData;

	int PropagateRank;
	double Mass;
	cplx TimeStep;

	void SetupLapackMatrices(const cplx &dt);
	void SetupBlasMatrices(const cplx &dt);

public:

	/*
	 * Get configuration options from configsection
	 */
	void ApplyConfigSection(const ConfigSection &config);

	/*
	 * Setup propagator object
	 */
	void Setup(const cplx &dt, const Wavefunction<Rank> &psi, BSpline::Ptr bsplineObject, int rank);

	/*
	 * Advance wavefunction one timestep
	 */
	void AdvanceStep(Wavefunction<Rank> &psi);

	/*
	 * Multiply hamilton matrix on wavefunction. Uses BLAS and LAPACK for speedy computation
	 */
	void MultiplyHamiltonian(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi);
	void ApplyCrankNicolson(const blitz::Array<cplx, 2> &matrix, blitz::Array<cplx, 3> &data);

	
	// Functions to return various propagator matrices
	blitz::Array<cplx, 2> GetOverlapMatrix() { return OverlapMatrix; }
	blitz::Array<cplx, 2> GetHamiltonianMatrix() { return HamiltonianMatrix; }
	blitz::Array<cplx, 2> GetPropagationMatrixBlas() { return PropagationMatrixBlas; }
	blitz::Array<cplx, 2> GetPropagationMatrix() { return PropagationMatrix; }
};

}; //Namespace

#endif


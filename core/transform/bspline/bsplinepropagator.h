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

	blitz::Array<cplx, 2> PropagationMatrix;
	blitz::Array<cplx, 2> HamiltonianMatrix;
	blitz::Array<cplx, 2> Eigenvectors;
	blitz::Array<double, 1> Eigenvalues;
	blitz::Array<cplx, 1> TempData;

	int PropagateRank;
	double Mass;

	void SetupHamiltonianMatrix(blitz::Array<cplx, 2> HamiltonianMatrix);
	void ComputeHamiltonianEigenvectors(blitz::Array<cplx, 2> HamiltonianMatrix);

public:

	void ApplyConfigSection(const ConfigSection &config);
	void Setup(const cplx &dt, const Wavefunction<Rank> &psi, int rank);
	void AdvanceStep(Wavefunction<Rank> &psi);
	void MultiplyHamiltonian(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi);
	void ApplyMatrix(const blitz::Array<cplx, 2> &matrix, blitz::Array<cplx, 3> &data);
	blitz::Array<cplx, 2> GetPropagationMatrix() { return PropagationMatrix; }
	blitz::Array<cplx, 2> GetHamiltonianMatrix() { return HamiltonianMatrix; }
	blitz::Array<cplx, 2> GetEigenvectors() { return Eigenvectors; }
	blitz::Array<double, 1> GetEigenvalues() { return Eigenvalues; }


};

}; //Namespace

#endif


#ifndef KRYLOVCOMMON_H
#define KRYLOVCOMMON_H

#include <core/common.h>
#include <core/utility/boostpythonhack.h>
#include "krylovbase.h"

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
	KrylovBase<Rank> *propagator = static_cast<KrylovBase<Rank>*>(data);
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
	if (real(outNorm) < 1e-10)
	{
		cout << "Outnorm is very small. Something is amiss." << endl;
		cout << "  InNorm = " << inNorm << endl;
		cout << " outNorm = " << outNorm << endl;
		cout << "  inAddr = " << inBuffer << endl;
		cout << " outAddr = " << outBuffer << endl;
		cout << "    diff = " << (int)(inBuffer - outBuffer) << endl;
	}
}



};

#endif


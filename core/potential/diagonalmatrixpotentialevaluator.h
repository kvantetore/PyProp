#ifndef DIAGONALMATRIXPOTENTIALEVALUATOR_H
#define DIAGONALMATRIXPOTENTIALEVALUATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../utility/blitzblas.h"

class DiagonalMatrixPotentialEvaluator
{
public:
	typedef Wavefunction<1>::DataArray VectorArray;

private:
	VectorArray DiagonalElement;
	VectorArray TempData;

public:
	
	void SetMatrixData(VectorArray diagonalElement)
	{
		this->DiagonalElement.resize(diagonalElement.shape());
		this->DiagonalElement = diagonalElement;
		this->TempData.resize(diagonalElement.extent(0));
	}
	
	void MultiplyPotential(Wavefunction<1> &psi, Wavefunction<1> &destPsi, double scaling)
	{
		VectorArray in = psi.GetData();
		VectorArray out = destPsi.GetData();

		out += in * scaling * DiagonalElement;
	}

	cplx CalculateExpectationValue(Wavefunction<1> &psi)
	{
		VectorArray in = psi.GetData();
		return sum(conj(in) * DiagonalElement * in);
	}
};


#endif


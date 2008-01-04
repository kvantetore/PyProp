#ifndef DENSEMATRIXPOTENTIALEVALUATOR_H
#define DENSEMATRIXPOTENTIALEVALUATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../utility/blitzblas.h"

class DenseMatrixPotentialEvaluator
{
public:
	typedef blitz::Array<cplx, 2> DataArray;
	typedef Wavefunction<1>::DataArray VectorArray;

private:
	DataArray MatrixElement;
	VectorArray TempData;

public:
	
	void SetMatrixData(DataArray matrixElement)
	{
		this->MatrixElement.resize(matrixElement.shape());
		this->MatrixElement = matrixElement;
		this->TempData.resize(matrixElement.extent(0));
	}
	
	void MultiplyPotential(Wavefunction<1> &psi, Wavefunction<1> &destPsi, double scaling)
	{
		VectorArray in = psi.GetData();
		VectorArray out = destPsi.GetData();
		
		TempData = 0;

		MatrixVectorMultiply(MatrixElement, in, TempData);

		out += scaling*TempData;

	}

	cplx CalculateExpectationValue(Wavefunction<1> &psi)
	{
		VectorArray in = psi.GetData();
		MatrixVectorMultiply(MatrixElement, in, TempData);
		return VectorInnerProduct(in, TempData);
	}
};


#endif


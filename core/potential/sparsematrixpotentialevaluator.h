#ifndef SPARSEMATRIXPOTENTIALEVALUATOR_H
#define SPARSEMATRIXPOTENTIALEVALUATOR_H

#include "../common.h"
#include "../wavefunction.h"


class SparseMatrixPotentialEvaluator
{
public:
	typedef blitz::Array<int, 1> IndexArray;
	typedef blitz::Array<cplx, 1> DataArray;

private:
	IndexArray RowIndex;
	IndexArray ColIndex;
	DataArray MatrixElement;

public:
	
	void SetMatrixData(IndexArray rowIndex, IndexArray colIndex, DataArray matrixElement)
	{
		this->RowIndex.resize(rowIndex.shape());
		this->RowIndex = rowIndex;
		this->ColIndex.resize(colIndex.shape());
		this->ColIndex = colIndex;
		this->MatrixElement.resize(matrixElement.shape());
		this->MatrixElement = matrixElement;
	}
	
	void MultiplyPotential(Wavefunction<1> &psi, Wavefunction<1> &destPsi, double scaling)
	{
		DataArray in(psi.GetData());
		DataArray out(destPsi.GetData());

		for (int i=0; i<RowIndex.size(); i++)
		{
			int r = RowIndex(i);
			int c = ColIndex(i);
			out(r) += scaling * MatrixElement(i) * in(c);
			out(c) += scaling * conj(MatrixElement(i)) * in(r);
		}

		//cout << "In = " << in << endl;
	}

	cplx CalculateExpectationValue(Wavefunction<1> &psi)
	{
		DataArray in = psi.GetData();
		cplx out = 0;

		for (int i=0; i<RowIndex.size(); i++)
		{
			int r = RowIndex(i);
			int c = ColIndex(i);
			out += conj(in(r)) * MatrixElement(i) * in(c);
			out += conj(in(c)) * conj(MatrixElement(i)) * in(r);
		}

		return out;
	}
};


#endif


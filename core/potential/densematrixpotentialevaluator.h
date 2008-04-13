#ifndef DENSEMATRIXPOTENTIALEVALUATOR_H
#define DENSEMATRIXPOTENTIALEVALUATOR_H

#include "../common.h"
#include "../wavefunction.h"
#include "../utility/blitzblas.h"
#include "../utility/blitzlapack.h"
#include "../utility/blitztricks.h"

/*
 * DenseMatrixPotentialEvaluator is a potential evaluator that is diagonal 
 * in all ranks but one. The rank which is not diagonal is MatrixRowRank.
 * The matrix elements may depend parametrically on the coordinate in 
 * all other ranks. 
 *
 * This can perhaps be thought of as the Born-Openheimer approximation, as 
 * we have decoupled one rank from the rest.
 *
 * If we have a Rank-dimensional wavefunction, we will then have a 
 * Rank+1-dimensional matrix. The current implementation requires
 * the column-rank of the matrix to be the last column. That way, all
 * indices is the same between the matrix and the wavefunction except 
 * for the last index in the matrix.
 *
 * If this is seen as a regular matrix, it can be seen as block diagonal
 * 
 * V = V(r_i; r_0, r_1, ...) = V(r_i(0)) X V(r_i(1)) X ...
 *
 * We can therefore diagonalize each block by itself
 *
 * exp(V) = exp(V(r_i(0))) X exp(V(r_i(1))) X ...
 *
 */

template<int Rank>
class DenseMatrixPotentialEvaluator
{
public:
	typedef blitz::Array<cplx, Rank+1> MatrixArray;
	typedef typename Wavefunction<Rank>::DataArray VectorArray;

private:
	MatrixArray MatrixElement;
	MatrixArray Eigenvectors;
	VectorArray Eigenvalues;
	int MatrixMultiplyAlgorithm;

	VectorArray TempData;

	int MatrixRowRank;		//The rank in the wavefunction which should be transformed
	int MatrixColRank;		//The rank in MatrixElement which is the columns in the matrix-vector product

public:

	DenseMatrixPotentialEvaluator() : MatrixMultiplyAlgorithm(0) {}

	MatrixArray GetEigenvectors()
	{
		return Eigenvectors;
	}

	VectorArray GetEigenvalues()
	{
		return Eigenvalues;
	}

	MatrixArray GetMatrixData()
	{
		return MatrixElement;
	}

	int GetAlgorithm()
	{
		return MatrixMultiplyAlgorithm;
	}

	void SetAlgorithm(int algo)
	{
		MatrixMultiplyAlgorithm = algo;
	}

	
	void SetMatrixData(MatrixArray matrixElement, int rowRank, int colRank)
	{
		if (rowRank >= colRank)
		{
			throw std::runtime_error("Error in SetMatrixData. rowRank should be less than colRank");
		}
		if (colRank != Rank)
		{
			throw std::runtime_error("Error in SetMatrixData. It is currently implemented only when colRank is the last rank");
		}
		if (rowRank != Rank-1)
		{
			cout << "WARNING: it is recommended that rowRank is the last rank of the problem." << endl;
		}
		if (matrixElement.extent(rowRank) != matrixElement.extent(colRank))
		{
			throw std::runtime_error("Error in SetMatrixData. Matrix must be square");
		}

		this->MatrixElement.resize(matrixElement.shape());
		this->MatrixElement = matrixElement;

		blitz::TinyVector<int, Rank> vectorShape;
		for (int i=0; i<Rank; i++)
		{
			vectorShape(i) = MatrixElement.extent(i);
		}
		this->TempData.resize(vectorShape);

		MatrixRowRank = rowRank;
		MatrixColRank = colRank;

		DiagonalizeMatrix();

		cout << matrixElement(0,5) << endl;
	}

	void DiagonalizeMatrix()
	{
		using namespace blitz;

		TinyVector<int, 2> shape;
		TinyVector<int, 2> stride;

		shape(0) = MatrixElement.extent(MatrixRowRank);
		shape(1) = MatrixElement.extent(MatrixColRank);
		stride(0) = MatrixElement.stride(MatrixRowRank);
		stride(1) = MatrixElement.stride(MatrixColRank);

		TinyVector<int, Rank> vectorIndex = 0;
		TinyVector<int, Rank> vectorShape;
		TinyVector<int, Rank+1> matrixIndex = 0;
		for (int i=0; i<Rank; i++)
		{
			vectorShape(i) = MatrixElement.extent(i);
		}

		//cout << "Matrix shape = " << MatrixElement.shape() << endl;
		//cout << "Vector shape = " << vectorShape << endl;

		Eigenvectors.resize(MatrixElement.shape());
		Eigenvalues.resize(vectorShape);

		Array<double, 1> eigenvalues(shape(0));
		linalg::LAPACK<cplx> lapack;

		//Iterate over all dimensions except MatrixRowRank and MatrixColRank
		int skipDimension = MatrixRowRank;   
		do 
		{
			//cout << "curIndex = " << ToString(vectorIndex) << endl;

			for (int i=0; i<Rank; i++)
			{
				matrixIndex(i) = vectorIndex(i);
			}

			//CalculateEigenvectorFactorization destroys the input matrix, so we make a copy
			//The matrix is in col-major storage, while CalculateEigenvectorFactorization expects 
			//row-major storage => Calculate left eigenvectors instead
			//A X = X L => X^T A^T = L X^T
			//(OTOH, the matrix should be Hermitian, so it doesn't matter much, as long as we get
			//the eigenvectors out correctly)
			Array<cplx, 2> matrixView(&MatrixElement(matrixIndex), shape, stride, neverDeleteData);
			Array<cplx, 2> curMatrix(matrixView.copy());
			lapack.CalculateEigenvectorFactorizationHermitian(true, linalg::LAPACK<cplx>::HermitianUpper, curMatrix, eigenvalues);

			//Store eigenvalues
			TinyVector<int, 1> eigenvalueShape(shape(0));
			TinyVector<int, 1> eigenvalueStride(Eigenvalues.stride(MatrixRowRank));
			Array<cplx, 1> eigenvalueView(&Eigenvalues(vectorIndex), eigenvalueShape, eigenvalueStride, neverDeleteData);
			eigenvalueView = eigenvalues(tensor::i);

			//Store eigenvectors
			Array<cplx, 2> eigenvectorView(&Eigenvectors(matrixIndex), shape, stride, neverDeleteData);
			eigenvectorView = curMatrix;

		} while (IncIndex(vectorShape, vectorIndex, Rank-1, skipDimension));

	}
	
	void MultiplyPotential(Wavefunction<Rank> &psi, Wavefunction<Rank> &destPsi, double scaling)
	{
		VectorArray in = psi.GetData();
		VectorArray out = destPsi.GetData();

		MultiplyPotential(in, out, scaling);
	}

	void MultiplyPotential(VectorArray &in, VectorArray &out, double scaling)
	{
		using namespace blitz;
		typedef typename MatrixArray::iterator iterator;
		
		TinyVector<int, Rank> inPos;
		TinyVector<int, Rank> outPos;
		for (iterator matrixIterator=MatrixElement.begin(); matrixIterator!=MatrixElement.end(); matrixIterator++)
		{
			for (int i=0; i<Rank; i++)
			{
				inPos(i) = matrixIterator.position()(i);
				outPos(i) = matrixIterator.position()(i);
			}
			inPos(MatrixRowRank) = matrixIterator.position()(Rank);
			
			out(outPos) += scaling * (*matrixIterator) * in(inPos);
		}

		//out += scaling * sum( MatrixElement(i, j, k) * in(i, k), k);
	}

	void ApplyPotential(Wavefunction<Rank> &psi, cplx dt, double scaling)
	{
		using namespace blitz;
		typedef typename MatrixArray::iterator iterator;

		VectorArray v = psi.GetData();

		//cout << "V = " << v << endl;

		//We have A = X L X*
		//=> exp(dt f(t) A) v = X exp(dt f(t) L) X* v
	
		//1) Transform into eigenspace of A
		if (MatrixMultiplyAlgorithm == 0)
		{
			TinyVector<int, Rank> inPos;
			TinyVector<int, Rank> outPos;
			TempData = 0;
			for (iterator matrixIterator=Eigenvectors.begin(); matrixIterator!=Eigenvectors.end(); matrixIterator++)
			{
				for (int i=0; i<Rank; i++)
				{
					inPos(i) = matrixIterator.position()(i);
					outPos(i) = matrixIterator.position()(i);
				}
				inPos(MatrixRowRank) = matrixIterator.position()(Rank);
				
				TempData(outPos) += *matrixIterator * v(inPos);
			}
		}
		if (MatrixMultiplyAlgorithm == 1)
		{
			TempData = sum(Eigenvectors(tensor::i, tensor::j, tensor::k) * v(tensor::i, tensor::k), tensor::k);
		}
		if (MatrixMultiplyAlgorithm == 2)
		{
			int N0 = v.extent(0);
			int N1 = v.extent(1);
			for (int i=0; i<N0; i++)
			{
				for (int j=0; j<N1; j++)
				{
					cplx value=0;
					for (int k=0; k<N1; k++)
					{
						value += Eigenvectors(i, j, k) * v(i, k);
					}
					TempData(i, j) = value;
				}
			}
		}

		//2) Apply temp := exp(df f(t) L) temp
		TempData *= exp(- cplx(0, 1) * Eigenvalues * dt * scaling);

		//3) Transform back to original space
		if (MatrixMultiplyAlgorithm == 0)
		{
			TinyVector<int, Rank> inPos;
			TinyVector<int, Rank> outPos;
			v = 0;
			for (iterator matrixIterator=Eigenvectors.begin(); matrixIterator!=Eigenvectors.end(); matrixIterator++)
			{
				for (int i=0; i<Rank; i++)
				{
					inPos(i) = matrixIterator.position()(i);
					outPos(i) = matrixIterator.position()(i);
				}
				outPos(MatrixRowRank) = matrixIterator.position()(Rank);
				
				v(outPos) += *matrixIterator * TempData(inPos);
			}
		}

		if (MatrixMultiplyAlgorithm == 1)
		{
			v = sum( Eigenvectors(tensor::i, tensor::k, tensor::j) * TempData(tensor::i, tensor::k), tensor::k);
		}

		if (MatrixMultiplyAlgorithm == 2)
		{
			int N0 = v.extent(0);
			int N1 = v.extent(1);
			for (int i=0; i<N0; i++)
			{
				for (int j=0; j<N1; j++)
				{
					cplx value=0;
					for (int k=0; k<N1; k++)
					{
						value += Eigenvectors(i, k, j) * TempData(i, k);
					}
					v(i, j) = value;
				}
			}
		}
		//cout << "v = " << v << endl;
	}

};


class DenseMatrixPotentialEvaluatorOld
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


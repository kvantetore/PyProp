#include "orthogonalpolynomial.h"
#include "fortran.h"

extern "C"
{
	void FORTRAN_NAME(dstev)(char* JOB, int *N, double* D, double *E, double *Z, int *LDZ, double *WORK, int *INFO);
};

using namespace blitz;

namespace OrthogonalPolynomial
{

Array<double, 2> CalculateQuadrature(int N, Array<double, 1> diagonal, Array<double, 1> subDiagonal)
{
	Array<double, 2> eigenVectors(N, N);
	Array<double, 1> work(2*N);
	char job = 'V';
	int info = 0;

	//Call diagonalization routine
	FORTRAN_NAME(dstev)(&job, &N, diagonal.data(), subDiagonal.data(), eigenVectors.data(), &N, work.data(), &info);

	//the zeros are the eigenvalues, and the weights are 2 * ((v_j)_0)**2 where v_j is the j'th eigenvector.
	//and why? because it works!
	Array<double, 2> quadrature(N, 2);
	quadrature(Range::all(), 0) = diagonal;
	quadrature(Range::all(), 1) = 2 * sqr(eigenVectors(Range::all(), 0));
	return quadrature;
}

Array<double, 2> CalculateQuadrature(int N, PolynomialType type)
{
	Array<double, 1> diagonal(N);
	Array<double, 1> subDiagonal(N-1);

	if (type == Legendre1)
	{
		diagonal = 0;
		Array<double, 1> temp(N-1);
		temp = 1.0 / ((tensor::i + 1) * (tensor::i + 1)); 
		subDiagonal = 1.0/sqrt(4. - temp);
	}
	else
	{
		throw std::runtime_error("Invalid polynomial type. Only legendre is currently supported");
	}

	return CalculateQuadrature(N, diagonal, subDiagonal);
}

} // namespace


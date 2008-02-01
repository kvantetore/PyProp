#ifndef PIRAM_FUNCTORS
#define PIRAM_FUNCTORS

#include <core/common.h>
#include <core/mpi/mpitraits.h>
#include "blitzblas.h"

namespace piram
{
using namespace blitz::linalg;

/* 
 * Abstract functor class to implement matrix-vector multiplication
 */
template<class T>
class OperatorFunctor
{
public:
	typedef boost::shared_ptr< OperatorFunctor > Ptr;

	virtual ~OperatorFunctor() {}
	virtual void operator()(blitz::Array<T, 1> &in, blitz::Array<T, 1> &out) = 0;
};

/*
 * Abstract functor class to implement initialization of residual
 */
template<class T>
class SetupResidualFunctor
{
public:
	typedef boost::shared_ptr< SetupResidualFunctor > Ptr;

	virtual ~SetupResidualFunctor() {}
	virtual void operator()(blitz::Array<T, 1> &residual) = 0;
};

/* Abstract functor class to implement calculation of norm and inner product
 *
 * Ok, so its not really a functor, as it doesnt overload operator(), but
 * if you please bear with me for a moment...
 */
template<class T, class NormType>
class IntegrationFunctor
{
public:
	typedef boost::shared_ptr< IntegrationFunctor > Ptr;
	typedef blitz::Array<T, 1> VectorType;
	typedef blitz::Array<T, 2> MatrixType;

	virtual ~IntegrationFunctor() {}

	virtual NormType Norm(VectorType &vector) = 0;
	virtual T InnerProduct(VectorType &left, VectorType &right) = 0;
	virtual void InnerProduct(MatrixType &left, VectorType &right, VectorType &out, VectorType &temp) = 0;
};

template<class T, class NormType>
class BLASIntegrationFunctor : public IntegrationFunctor<T, NormType>
{
public:
	typedef boost::shared_ptr< BLASIntegrationFunctor > Ptr;
	typedef blitz::Array<T, 1> VectorType;
	typedef blitz::Array<T, 2> MatrixType;

private:
	typedef blitz::linalg::BLAS<T> BLAS;

	bool DisableMPI;
	MPI_Comm CommBase;
	BLAS blas;

	void GlobalAddVector(VectorType &in, VectorType &out)
	{
		if (DisableMPI)
		{
			out = in;
		}

		MPITraits<T> traits;
		MPI_Allreduce(in.data(), out.data(), in.extent(0)*traits.Length(), traits.Type(), MPI_SUM, CommBase);
	}

public:
	BLASIntegrationFunctor(bool disableMPI, MPI_Comm commBase) : DisableMPI(disableMPI), CommBase(commBase) {}
	virtual ~BLASIntegrationFunctor() {}

	virtual NormType Norm(VectorType &vector)
	{
		return sqrt((NormType) InnerProduct(vector, vector).real());
	}

	virtual T InnerProduct(VectorType &left, VectorType &right) 
	{
		T localValue = blas.InnerProduct(left, right);
		if (DisableMPI)
		{
			return localValue;
		}

		T globalValue;
		MPITraits<T> traits;
		MPI_Allreduce(&localValue, &globalValue, 1*traits.Length(), traits.Type(), MPI_SUM, CommBase);

		return globalValue;
	}

	virtual void InnerProduct(MatrixType &left, VectorType &right, VectorType &out, VectorType &temp)
	{
		blas.MultiplyMatrixVector(MatrixTranspose::Conjugate, left, 1.0, right, 0.0, temp);
		GlobalAddVector(temp, out);
	}
};


} //Namespace

#endif


#ifndef PIRAM_FUNCTORS
#define PIRAM_FUNCTORS

#include <core/common.h>
#include <core/mpi/mpitraits.h>
#include "blitzblas.h"

#include <functional>
#include <algorithm>

namespace piram
{
using namespace blitz::linalg;

/* 
 * Abstract functor class to decide which shifts to use in the 
 * implicit restarting process
 */
template<class T>
class ShiftFunctor
{
public:
	typedef boost::shared_ptr< ShiftFunctor > Ptr;

	virtual ~ShiftFunctor() {}
	virtual void operator()(blitz::Array<T, 1> eigenvalues, blitz::Array<int, 1> ordering) = 0;
};

/*
 * Implementation of the ShiftFunctor for using the least desirable ritz values
 * as shifts. 
 */
template<class T>
class CompareType : public std::binary_function<T, T, bool>
{
public:
	typedef shared_ptr< CompareType<T> > Ptr;

	virtual bool operator()(T a, T b) = 0;
};

template<class T, class CompareT>
class ShiftFunctorSelectRitzValues : public ShiftFunctor<T>
{
public:
	ShiftFunctorSelectRitzValues(const CompareT &compare) : Compare(compare) {}
	
	virtual void operator()(blitz::Array<T, 1> ritzvalues, blitz::Array<int, 1> ordering) 
	{
		/*
		//Sort ritz values such that the shifts are in the last end of the array
		//cout << "E = " << real(ritzvalues) << endl;
		
		std::sort(ritzvalues.data(), ritzvalues.data() + ritzvalues.extent(0), Compare);
		//cout << "sort(E) = " << real(ritzvalues) << endl;

		int startIndex = ritzvalues.extent(0) - shifts.extent(0);
		shifts = ritzvalues(blitz::Range(startIndex, blitz::Range::toEnd));
		//cout << "shifts(E) = " << real(shifts) << endl;
		//cout << endl;
		*/

		//TODO: For f*** sake, don't implement our own sort algo.  
		
		int n = ritzvalues.extent(0);
		for (int i=0; i<n; i++)
		{
			for (int j=i+1; j<n; j++)
			{
				if (! Compare( ritzvalues(i), ritzvalues(j)) )
				{
					T temp = ritzvalues(i);
					ritzvalues(i) = ritzvalues(j);
					ritzvalues(j) = temp;
		
					int itemp = ordering(i);
					ordering(i) = ordering(j);
					ordering(j) = itemp;
				}
			}
		}

	}

private:
	CompareT Compare;
};


/* 
 *   Compare functor for the ShiftFunctorSelectRitzValues which chooses the 
 *   highest real valued ritz values as shifts. This makes the ritz values
 *   converge towards the lowest eigenvalues of the system
 */
class CompareComplexLessReal : public CompareType<cplx>
{
public:
	virtual bool operator()(cplx a, cplx b)
	{
		return real(a) < real(b);
	}
};

/*
 * Compare functor for converging towards the eigenvalues close to Shift
 */
class CompareComplexNearShift : public CompareType<cplx>
{
public:
	CompareComplexNearShift(cplx shift) : Shift(shift) { }
	virtual bool operator()(cplx a, cplx b)
	{
		return std::abs(a-Shift) < std::abs(b-Shift);
	}

private:
	cplx Shift; 
};


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
		else
		{
			MPITraits<T> traits;
			MPI_Allreduce(in.data(), out.data(), in.extent(0)*traits.Length(), traits.Type(), MPI_SUM, CommBase);
		}
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
		else
		{
			T globalValue;
			MPITraits<T> traits;
			MPI_Allreduce(&localValue, &globalValue, 1*traits.Length(), traits.Type(), MPI_SUM, CommBase);

			return globalValue;
		}
	}

	virtual void InnerProduct(MatrixType &left, VectorType &right, VectorType &out, VectorType &temp)
	{
		blas.MultiplyMatrixVector(MatrixTranspose::Conjugate, left, 1.0, right, 0.0, temp);
		GlobalAddVector(temp, out);
	}
};


} //Namespace

#endif


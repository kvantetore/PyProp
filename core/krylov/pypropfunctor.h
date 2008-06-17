#include <core/common.h>
#include <core/wavefunction.h>
#include <core/mpi/distributedmodel.h>
#include <core/utility/blitztricks.h>

class PypropKrylovWrapper
{
public:
	virtual ~PypropKrylovWrapper() {}

	void virtual SetupResidual(blitz::Array<cplx, 1> &residual);
	void virtual ApplyOperator(blitz::Array<cplx, 1> &input, blitz::Array<cplx, 1> &output);
};

using namespace boost::python;

/* Functors for pIRAM */
template<int Rank>
class PypropOperatorFunctor : public piram::OperatorFunctor<cplx>
{
public:
	typedef boost::shared_ptr< PypropOperatorFunctor > Ptr;

	PypropOperatorFunctor(PypropKrylovWrapper *solver) : Solver(solver) {}

	virtual void operator()(blitz::Array<cplx, 1> &in, blitz::Array<cplx, 1> &out) 
	{
		Solver->ApplyOperator(in, out);
	}

private:
	PypropKrylovWrapper* Solver;
};

template<int Rank>
class PypropSetupResidualFunctor : public piram::SetupResidualFunctor<cplx>
{
public:
	typedef boost::shared_ptr< PypropSetupResidualFunctor > Ptr;

	PypropSetupResidualFunctor(PypropKrylovWrapper *solver) : Solver(solver) {}

	virtual void operator()(blitz::Array<cplx, 1> &residual) 
	{
		Solver->SetupResidual(residual);
	}

private:
	PypropKrylovWrapper* Solver;
};

template<int Rank>
class PypropIntegrationFunctor : public piram::IntegrationFunctor<cplx, double>
{
public:
	typedef boost::shared_ptr< PypropIntegrationFunctor<Rank> > Ptr;
	typedef blitz::Array<cplx, 1> VectorType;
	typedef blitz::Array<cplx, 2> MatrixType;
	typedef typename Wavefunction<Rank>::DataArray DataArray;

private:
	blitz::TinyVector<int, Rank> Shape;
	blitz::TinyVector<int, Rank> Stride;
	typename Wavefunction<Rank>::Ptr Psi1;
	typename Wavefunction<Rank>::Ptr Psi2;


public:
	PypropIntegrationFunctor(typename Wavefunction<Rank>::Ptr psi1, typename Wavefunction<Rank>::Ptr psi2) : Psi1(psi1), Psi2(psi2) 
	{
		Shape = Psi1->GetData().shape();
		Stride = Psi1->GetData().stride();
	}

	virtual ~PypropIntegrationFunctor() {}

	virtual double Norm(VectorType &vector)
	{
		//Map the 1d vector to a a blitz array of correct shape
		DataArray array(vector.data(), Shape, Stride, blitz::neverDeleteData);

		//Set psi to point to correct buffer
		DataArray oldArray(Psi1->GetData());
		Psi1->SetData(array);

		double norm = Psi1->GetNorm();

		//restore the old buffer
		Psi1->SetData(oldArray);

		return norm;
	}

	virtual cplx InnerProduct(VectorType &left, VectorType &right) 
	{
		//Map the 1d vector to a a blitz array of correct shape
		DataArray leftArray(left.data(), Shape, Stride, blitz::neverDeleteData);
		DataArray rightArray(right.data(), Shape, Stride, blitz::neverDeleteData);

		//Set psi to point to correct buffer
		DataArray oldArray1(Psi1->GetData());
		DataArray oldArray2(Psi2->GetData());
		Psi1->SetData(leftArray);
		Psi2->SetData(rightArray);

		cplx value = Psi2->InnerProduct(*Psi1);

		//Restore the old buffers
		Psi1->SetData(oldArray1);
		Psi2->SetData(oldArray2);

		return value;
	}

	virtual void InnerProduct(MatrixType &left, VectorType &right, VectorType &out, VectorType &temp)
	{
		DataArray oldArray1(Psi1->GetData());
		DataArray oldArray2(Psi2->GetData());
		
		//Set psi to point to correct buffer
		DataArray rightArray(right.data(), Shape, Stride, blitz::neverDeleteData);
		Psi2->SetData(rightArray);

		//Perform local inner products
		int vectorCount = left.extent(0);
		for (int i=0; i<vectorCount; i++)
		{
			DataArray leftArray(&left(i,0), Shape, Stride, blitz::neverDeleteData);
			Psi1->SetData(leftArray);

			temp(i) = Psi2->LocalInnerProduct(*Psi1);
		}

		//Do the sum over all processors
		Psi1->GetRepresentation()->GetDistributedModel()->GetGlobalSum(temp, out);

		//Restore the old buffers
		Psi1->SetData(oldArray1);
		Psi2->SetData(oldArray2);
	}
};




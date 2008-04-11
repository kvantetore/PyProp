#include "pampwrapper.h"

#include <core/wavefunction.h>
#include <core/representation/representation.h>
#include <core/mpi/distributedmodel.h>
#include <core/utility/fortran.h>
#include <core/utility/blitztricks.h>

#include "../piram/piram/functors.h"

namespace krylov
{

using namespace boost::python;

/* Functors for pAMP */
template<int Rank>
class PypropOperatorFunctor : public piram::OperatorFunctor<cplx>
{
public:
	typedef boost::shared_ptr< PypropOperatorFunctor > Ptr;

	PypropOperatorFunctor(PampWrapper<Rank> *propagator) : Propagator(propagator) {}

	virtual void operator()(blitz::Array<cplx, 1> &in, blitz::Array<cplx, 1> &out) 
	{
		Propagator->ApplyOperator(in, out);
	}

private:
	PampWrapper<Rank>* Propagator;
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


/* Implementation of PampWrapper */
template<int Rank>
void PampWrapper<Rank>::ApplyOperator(blitz::Array<cplx, 1> &input, blitz::Array<cplx, 1> &output)
{
	if (Psi == 0)
	{
		throw std::runtime_error("Psi is 0");
	}
	if (TempPsi == 0)
	{
		throw std::runtime_error("TempPsi is 0");
	}

	//Map the 1d vectors to a a blitz array of correct shape
	DataVector shape = Psi->GetData().shape();
	DataVector stride = Psi->GetData().stride();
	DataArray inData(input.data(), shape, stride, blitz::neverDeleteData);
	DataArray outData(output.data(), shape, stride, blitz::neverDeleteData);
	outData = 0;

	//Set psi and tempPsi to point to correct vectorsbuffers
	DataArray oldData = Psi->GetData();
	DataArray oldTempData = TempPsi->GetData();
	Psi->SetData(inData);
	TempPsi->SetData(outData);

	Callback(Psi, TempPsi, TimeStep, CurTime);

	//Restore the former buffers
	Psi->SetData(oldData);
	TempPsi->SetData(oldTempData);
}


template<int Rank>
void PampWrapper<Rank>::ApplyConfigSection(const ConfigSection &config)
{
	config.Get("krylov_basis_size", Propagator.BasisSize);
}


template<int Rank>
void PampWrapper<Rank>::Setup(const typename Wavefunction<Rank>::Ptr psi)
{
	Propagator.MatrixSize = psi->GetData().size();
    Propagator.DisableMPI = psi->GetRepresentation()->GetDistributedModel()->IsSingleProc();
	Propagator.MatrixOperator = typename PypropOperatorFunctor<Rank>::Ptr( new PypropOperatorFunctor<Rank>(this) );
	Propagator.Setup();
}


template<int Rank>
void PampWrapper<Rank>::AdvanceStep(object callback, typename Wavefunction<Rank>::Ptr psi, typename Wavefunction<Rank>::Ptr tempPsi, cplx dt, double t)
{
	//The callback-functions uses these variables
	this->Psi = psi;
	this->TempPsi = tempPsi;
	this->Callback = callback;
	this->CurTime = t;
	this->TimeStep = dt;

	//Use our custom integration
	typename piram::IntegrationFunctor<cplx, double>::Ptr integration(new PypropIntegrationFunctor<Rank>(Psi, TempPsi));
	Propagator.Integration = integration;

	typename Wavefunction<Rank>::DataArray vector = psi->GetData();
	blitz::Array<cplx, 1> vector1d = MapToRank1(vector);
	Propagator.PropagateVector(vector1d, dt);

	//Zero the pointers to avoid mishaps
	this->Psi = typename Wavefunction<Rank>::Ptr();
	this->TempPsi = typename Wavefunction<Rank>::Ptr();
}

template class PampWrapper<1>;
template class PampWrapper<2>;
template class PampWrapper<3>;
template class PampWrapper<4>;

} //Namespace
	

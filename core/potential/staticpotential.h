#ifndef STATICPOTENTIAL_H
#define STATICPOTENTIAL_H

#include "../common.h"
#include "../wavefunction.h"
#include "../utility/blitzblas.h"

template<int Rank>
class StaticPotential
{
public:
	enum StorageModel
	{
		StorageValue = 1,
		StorageExpValue = 2, 
		StorageBoth = 1+2
	};

private:
	StorageModel Storage;
	blitz::Array<cplx, Rank> PotentialData;
	blitz::Array<cplx, Rank> PotentialDataExp;

public:
	StaticPotential() {}
	~StaticPotential() {}

	bool UseStorageValue()
	{
		return Storage == StorageValue || Storage == StorageBoth;
	}
	
	bool UseStorageExpValue()
	{
		return Storage == StorageExpValue || Storage == StorageBoth;
	}

	void InitializePotential(Wavefunction<Rank> &psi, StorageModel storage)
	{
		Storage = storage;
	
		int storageMultiply = 1;
		if (Storage == StorageBoth)
		{
			storageMultiply = 2;
		}
		std::cout << "Allocating StaticPotential of shape " << storageMultiply << " * " << psi.Data.shape() 
		          << " (~" << blitz::product(psi.Data.shape()) * storageMultiply * sizeof(cplx) / (1024*1024) <<
			  "MB)" 
			  << std::endl;

		if (UseStorageValue())
		{
			PotentialData.resize(psi.Data.shape());
		}
		if (UseStorageExpValue())
		{
			PotentialDataExp.resize(psi.Data.shape());
		}

	}
	
	blitz::Array<cplx, Rank> GetPotentialData()
	{
		return PotentialData;
	}

	blitz::Array<cplx, Rank> GetPotentialDataExp()
	{
		return PotentialDataExp;
	}

	StorageModel GetStorageModel()
	{
		return Storage;
	}
	
	void ApplyPotential(Wavefunction<Rank> &psi, cplx dt, const cplx scaling)
	{
		ValidatePsi(psi);

		if (UseStorageExpValue())
		{
			VectorElementMultiply(psi.Data, PotentialDataExp, psi.Data);
		}
		else
		{
			const cplx timeFactor = - cplx(0.0, 1.0) * dt;
			psi.GetData() *= exp(timeFactor * scaling * PotentialData);
		}
	}

	void MultiplyPotential(Wavefunction<Rank> &psi, Wavefunction<Rank> &destPsi, const cplx &dt, const cplx &scaling)
	{
		ValidatePsi(psi);

		typename Wavefunction<Rank>::DataArray dest(destPsi.GetData());
		typename Wavefunction<Rank>::DataArray src(psi.GetData());

		if (UseStorageValue())
		{
			dest += scaling * PotentialData * src;
		}
		else
		{
			const cplx imaginaryUnit = cplx(0.0, 1.0);
			cplx scale = - 1.0 / (imaginaryUnit * dt);
			dest += scaling * log(PotentialDataExp) * scale * src;
		}
	}


private:
	void ValidatePsi(const Wavefunction<Rank> &psi)
	{
		blitz::TinyVector<int, Rank> potentialShape;
		blitz::TinyVector<int, Rank> potentialOrdering;
		if (UseStorageValue()) 
		{
			potentialShape = PotentialData.shape();
		}
		else 
		{
			potentialShape = PotentialDataExp.shape();
		}

		blitz::TinyVector<int, Rank> psiShape = psi.Data.shape();
		if (potentialShape != psiShape)
		{
			std::cout << "Error: Potential has different shape than Wavefunction: " 
			          << PotentialData.shape() << " != " << psi.Data.shape() << std::endl;
			     
			exit(-1);
		}
	}
};

#endif


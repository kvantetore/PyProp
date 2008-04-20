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
		StorageExpValue = 2
	};

private:
	StorageModel Storage;
	blitz::Array<cplx, Rank> PotentialData;

public:
	StaticPotential() {}
	
	void InitializePotential(Wavefunction<Rank> &psi, StorageModel storage)
	{
		std::cout << "Allocating StaticPotential of shape " << psi.Data.shape() 
		          << " (~" << blitz::product(psi.Data.shape()) * sizeof(cplx) / (1024*1024) <<
			  "MB)" 
			  << std::endl;
		PotentialData.resize(psi.Data.shape());
		Storage = storage;
	}
	
	blitz::Array<cplx, Rank> GetPotentialData()
	{
		return PotentialData;
	}

	StorageModel GetStorageModel()
	{
		return Storage;
	}
	
	void ApplyPotential(Wavefunction<Rank> &psi, cplx dt, const cplx scaling)
	{
		ValidatePsi(psi);

		if (Storage == StorageExpValue)
		{
			VectorElementMultiply(psi.Data, PotentialData, psi.Data);
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

		if (Storage == StorageExpValue)
		{
			const cplx imaginaryUnit = cplx(0.0, 1.0);
			cplx scale = - 1.0 / (imaginaryUnit * dt);
			dest += scaling * log(PotentialData) * scale * src;
		}
		else
		{
			dest += scaling * PotentialData * src;
		}
	}


private:
	void ValidatePsi(const Wavefunction<Rank> &psi)
	{
		blitz::TinyVector<int, Rank> potentialShape = PotentialData.shape();
		blitz::TinyVector<int, Rank> psiShape = psi.Data.shape();
		if (potentialShape != psiShape)
		{
			std::cout << "Error: Potential has different shape than Wavefunction: " 
			          << PotentialData.shape() << " != " << psi.Data.shape() << std::endl;
			     
			exit(0);
		}
		if (psi.Data.ordering() != PotentialData.ordering())
		{
			std::cout << "Warning: Potential has different ordering than Wavefunction. "
				      << "Possible performance penalty. " << std::endl;
		}
	}
};

#endif


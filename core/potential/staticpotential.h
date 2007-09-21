#ifndef STATICPOTENTIAL_H
#define STATICPOTENTIAL_H

#include "../common.h"
#include "../wavefunction.h"
#include "../utility/blitzblas.h"

template<int Rank>
class StaticPotential
{
private:
	//Exponential form of the potential:
	//PotentialData(x) = exp(- i * dt * V(x) )
	blitz::Array<cplx, Rank> PotentialData;

public:

	StaticPotential() {}
	
	void InitializePotential(Wavefunction<Rank> &psi)
	{
		std::cout << "Allocating StaticPotential of shape " << psi.Data.shape() 
		          << " (~" << blitz::product(psi.Data.shape()) * sizeof(cplx) / (1024*1024) <<
			  "MB)" 
			  << std::endl;
		PotentialData.changeOrdering(psi.Data.ordering());
		PotentialData.resize(psi.Data.shape());
	}
	
	blitz::Array<cplx, Rank> GetPotentialData()
	{
		return PotentialData;
	}
	
	void ApplyPotential(Wavefunction<Rank> &psi)
	{
		ValidatePsi(psi);

		VectorElementMultiply(psi.Data, PotentialData, psi.Data);
	}

	void MultiplyPotential(Wavefunction<Rank> &psi, Wavefunction<Rank> &destPsi, const cplx &dt)
	{
		ValidatePsi(psi);
		const cplx imaginaryUnit = cplx(0.0, 1.0);
		cplx scale = - 1.0 / (imaginaryUnit * dt);
		typename Wavefunction<Rank>::DataArray dest(destPsi.GetData());
		typename Wavefunction<Rank>::DataArray src(psi.GetData());
		dest += log(PotentialData) * scale * src;
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

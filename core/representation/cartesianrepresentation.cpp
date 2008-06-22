#include "cartesianrepresentation.h"
#include "../wavefunction.h"
#include "../utility/blitzblas.h"

//Implementation of the Representation interface.
template<int Rank>
blitz::TinyVector<int, Rank> CartesianRepresentation<Rank>::GetFullShape()
{
	blitz::TinyVector<int, Rank> fullShape;
	for (int i=0;i<Rank;i++)
	{
		fullShape(i) = Range(i).Count;
	}
	return fullShape;
}

template<int Rank>
cplx CartesianRepresentation<Rank>::InnerProduct(const Wavefunction<Rank> &w1, const Wavefunction<Rank> &w2)
{
	double weight = GetScalarWeight();
	return VectorInnerProduct(w1.Data, w2.Data) * weight;
}

template <int Rank>
void CartesianRepresentation<Rank>::ApplyConfigSection(const ConfigSection &cfg)
{
	//Check that the rank specified in the config file is prop
	int configRank = 0;
	cfg.Get("rank", configRank);
	
	if (configRank != Rank)
	{
		std::cout << "Invalid rank (" << configRank << ") "
				<< "specified for CartesianRange<" << Rank << ">" << std::endl;
		exit(-1);
	}

	//Get range	
	for (int i=0; i<Rank; i++)
	{
		blitz::TinyVector<double, 3> rangeVector;
		std::string settingName = "rank" + ToString(i);
		cfg.Get(settingName, rangeVector);
		
		double min = rangeVector(0);
		double max = rangeVector(1);
		int count = (int)rangeVector(2);
		
		Range(i) = CartesianRange(min, max, count);
	}
}

template<int Rank>
void CartesianRepresentation<Rank>::MultiplyOverlap(Wavefunction<Rank> &srcPsi, Wavefunction<Rank> &dstPsi, int rank)
{
	int activeRank = rank - this->GetBaseRank();
	if (activeRank < 0 || activeRank >= Rank)
	{
		cout << "Rank " << rank << " invalid for CartesianRepresntation" << endl;
		throw std::runtime_error("Invalid Rank");
	}
	CopyVector(Range(activeRank).Dx, srcPsi.GetData(), 0.0, dstPsi.GetData());
}

template<int Rank>
void CartesianRepresentation<Rank>::MultiplyOverlap(Wavefunction<Rank> &psi)
{
	ScaleVector(GetScalarWeight(), psi.GetData());
}

template<int Rank>
void CartesianRepresentation<Rank>::SolveOverlap(Wavefunction<Rank> &psi)
{
	ScaleVector(1./GetScalarWeight(), psi.GetData());
}

template<int Rank>
void CartesianRepresentation<Rank>::MultiplySqrtOverlap(bool conjugate, Wavefunction<Rank> &psi)
{
	ScaleVector(sqrt(GetScalarWeight()), psi.GetData());
}

template<int Rank>
void CartesianRepresentation<Rank>::SolveSqrtOverlap(bool conjugate, Wavefunction<Rank> &psi)
{
	ScaleVector(1.0/sqrt(GetScalarWeight()), psi.GetData());
}

template<int Rank>
OverlapMatrix::Ptr CartesianRepresentation<Rank>::GetGlobalOverlapMatrix(int rank)
{
	int activeRank = rank - this->GetBaseRank();
	if (activeRank < 0 || activeRank >= Rank)
	{
		cout << "Rank " << rank << " invalid for CartesianRepresntation" << endl;
		throw std::runtime_error("Invalid Rank");
	}
	return Range(activeRank).GetOverlapMatrix();
}

template class CartesianRepresentation<1>;
template class CartesianRepresentation<2>;
template class CartesianRepresentation<3>;
template class CartesianRepresentation<4>;


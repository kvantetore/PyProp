#include <core/common.h>
#include <core/wavefunction.h>

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

using namespace boost::python;
using namespace blitz;

typedef Array<cplx, 2> MatrixType;
typedef Array<cplx, 1> VectorType;

list CalculatePopulationRadialProductStates(int l1, MatrixType V1, int l2, MatrixType V2, Array<cplx, 3> psiData, Array<int, 1> angularIndices)
{
#ifndef __GCCXML__ //GCCXML does not deal well with boost::python code
	list popList;

	int count0 = angularIndices.extent(0);
	int count1 = V1.extent(1);
	int count2 = V2.extent(1);
	int rcount = V1.extent(0);

	blitz::Array<cplx, 3> proj(count0, count1, count2);

	for (int i0=0; i0<angularIndices.extent(0); i0++)
	{
		int angIdx = angularIndices(i0);
		MatrixType psiSlice = psiData(angIdx, Range::all(), Range::all());

		for (int r1=0; r1<rcount; r1++)
		{
			for (int r2=0; r2<rcount; r2++)
			{
				for (int i1=0; i1<count1; i1++)
				{
					for (int i2=0; i2<count2; i2++)
					{
						cplx symFunc = V1(r1, i1) * V2(r2, i2);
						proj(i0, i1, i2) += conj(symFunc) * psiSlice(r1, r2);
					}
				}
			}
		}
	}

	proj *= conj(proj);

	for (int i1=0; i1<count1; i1++)
	{
		for (int i2=0; i2<count2; i2++)
		{
			double pop = real(sum(proj(Range::all(), i1, i2)));
			popList.append( make_tuple(i1, i2, pop));
		}
	}
	return popList;
#endif //__GCCXML__
}

Wavefunction<3>::Ptr GetWavefunctionParticleExchange(Wavefunction<3>::Ptr psi, list angularSymmetrizationPairs)
{
	typedef Array<cplx, 3> ArrayType;
	ArrayType data = psi->GetData();

	int countr = data.extent(1);
	typedef stl_input_iterator<tuple> Iterator;
	Iterator begin(angularSymmetrizationPairs);
	Iterator end;

	Wavefunction<3>::Ptr exchgPsi = psi->Copy();
	ArrayType exchgData = exchgPsi->GetData();

	for (Iterator i=begin; i!=end; i++)
	{
		int a1 = extract<int>((*i)[0]);
		int a2 = extract<int>((*i)[1]);

		for (int r1=0; r1<countr; r1++)
		{
			for (int r2=0; r2<countr; r2++)
			{
				exchgData(a1, r1, r2) = data(a2, r2, r1);
			}
		}
	}
	
	return exchgPsi;
}

void SymmetrizeWavefunction(Wavefunction<3>::Ptr psi, list angularSymmetrizationPairs, bool symmetrize)
{
	double symFactor = (symmetrize) ? 1 : -1;
	double normFactor = 0.5;

	Wavefunction<3>::Ptr exchgPsi = GetWavefunctionParticleExchange(psi, angularSymmetrizationPairs);
	psi->GetData() += symFactor * exchgPsi->GetData();
	psi->GetData() *= normFactor;
}

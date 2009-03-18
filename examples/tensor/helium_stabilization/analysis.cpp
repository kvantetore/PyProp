#include <core/common.h>
#include <core/utility/boostpythonhack.h>

using namespace blitz;

typedef Array<cplx, 2> MatrixType;
typedef Array<cplx, 1> VectorType;

list CalculatePopulationRadialProductStates(int l1, MatrixType V1, int l2, MatrixType V2, blitz::Array<cplx, 3> psiData, blitz::Array<int, 1> angularIndices)
{
#ifndef __GCCXML__ //GCCXML does not deal well with boost::python code
	list popList;

	int count0 = angularIndices.extent(0);
	int count1 = V1.extent(0);
	int count2 = V2.extent(0);
	int rcount = V1.extent(1);

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
						proj(i0, i1, i2) += conj(V1(r1, i1)) * conj(V2(r2, i2)) * psiSlice(r1, r2);
					}
				}
			}
		}
	}

	return popList;
#endif //__GCCXML__
}



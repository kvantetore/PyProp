#include "bsplinetransform.h"
#include "../../utility/blitztricks.h"

template<int Rank> class CombinedRepresentation;

namespace BSpline
{

template<int Rank>
void BSplineTransform<Rank>::SetupStep(Wavefunction<Rank> &psi, BSpline::Ptr bsplineObject, int baseRank)
{
	using namespace blitz;
	
	SetBaseRank(baseRank);
	BSplineObject = bsplineObject;

	/*
	 * Get shape of wavefunction (we are in the bspline repr). 
	 * Then allocate bspline grid representation buffer of wavefunction
	 * and store buffer names on object.
	 */
	TinyVector<int, Rank> gridShape = psi.GetData().shape();
	gridShape(baseRank) = BSplineObject->GetQuadratureGridGlobal().extent(0);
	BSplineGridDataName = psi.AllocateData(gridShape);
	BSplineDataName = psi.GetActiveBufferName();

	TempData.resize(BSplineObject->GetQuadratureGridGlobal().extent(0));
}

/*
 * Transform wavefunction from grid representation to the b-spline
 * basis representation. We call on a BSPLINE-class function to do
 * the actual 1D transforms. This involves solving a generalized
 * system of equations, using LAPACK.
 */
template<int Rank> 
void BSplineTransform<Rank>::ForwardTransform(Wavefunction<Rank> &psi)
{
	using namespace blitz;

	if(psi.GetActiveBufferName() != BSplineGridDataName)
	{
		throw std::runtime_error("Active databuffer is not what is should be...");
	}

	//Project input and output into 3D arrays with bsplineRank in the middle. 
	Array<cplx, Rank> srcData(psi.GetData());
	Array<cplx, Rank> dstData(psi.GetData(BSplineDataName));
	Array<cplx, 3> input3d = MapToRank3(srcData, BaseRank, 1);
	Array<cplx, 3> output3d = MapToRank3(dstData, BaseRank, 1);

	int psiSliceExtent = input3d.extent(1);
	Array<cplx, 1> psiSlice(psiSliceExtent);

	output3d = 0;
	int preCount = input3d.extent(0);
	int postCount = input3d.extent(2);
	for (int i = 0; i < preCount; i++)
	{
		for (int j = 0; j < postCount; j++)
		{
			// View slice of psi along propagation direction
			psiSlice = input3d(i, Range::all(), j);

			// Call on BSpline function to perform expansion
			TempData = 0;
			BSplineObject->ExpandFunctionInBSplines(psiSlice, TempData);
			output3d(i, Range::all() , j) = TempData;
		}
	}

	psi.SetActiveBuffer(BSplineDataName);
}

/*
 *
 */
template<int Rank>
void BSplineTransform<Rank>::InverseTransform(Wavefunction<Rank> &psi)
{
	using namespace blitz;

	if(psi.GetActiveBufferName() != BSplineDataName)
	{
		throw std::runtime_error("Active databuffer is not what is should be...");
	}

	//Project input and output into 3D arrays with bsplineRank in the middle. 
	Array<cplx, Rank> inputData(psi.GetData());
	Array<cplx, Rank> outputData(psi.GetData(BSplineGridDataName));
	Array<cplx, 3> input3d = MapToRank3(inputData, BaseRank, 1);
	Array<cplx, 3> output3d = MapToRank3(outputData, BaseRank, 1);

	int psiSliceExtent = input3d.extent(1);
	Array<cplx, 1> psiSlice(psiSliceExtent);

	output3d = 0;
	int preCount = input3d.extent(0);
	int postCount = input3d.extent(2);
	for (int i=0; i<preCount; i++)
	{
		for (int j=0; j<postCount; j++)
		{
			/*
			 * Copy wavefunction slice along bspline rank to temp
			 * array.
			 */
			psiSlice = input3d(i, Range::all(), j).copy();

			// Call on BSpline function to perform expansion
			output3d(i, Range::all() , j) = 
				BSplineObject->ConstructFunctionFromBSplineExpansion(psiSlice);
		}
	}
	psi.SetActiveBuffer(BSplineGridDataName);
}

template class BSplineTransform<1>;
template class BSplineTransform<2>;
template class BSplineTransform<3>;
template class BSplineTransform<4>;

} // Namespace

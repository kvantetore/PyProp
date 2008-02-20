#include "bsplinetransform.h"
#include "bspline.h"
#include "../../utility/blitztricks.h"

namespace BSpline
{

void BSplineTransform::SetupStep()
{

}

/*
 * Transform wavefunction from grid representation to the b-spline
 * basis representation. We call on a BSPLINE-class function to do
 * the actual 1D transforms. This involves solving a generalized
 * system of equations, using LAPACK.
 */
template<int Rank> 
void BSplineTransform::ForwardTransform(Wavefunction<Rank> &psi)
{
	using namespace blitz;

	BSpline::Ptr bsplineObject = 
		psi.GetRepresentation().GetRepresentation(BaseRank).GetBSplineObject();

	//Project input and output into 3D arrays with bsplineRank in the middle. 
	Array<cplx, 3> input3d = MapToRank3(psi, psi.Rank - 1, 1);
	Array<cplx, 3> output3d = MapToRank3(output, psi.Rank - 1, 1);

	int psiSliceExtent = input3d.extent(1);
	Array<cplx, 1> psiSlice(psiSliceExtent);

	output3d = 0;
	int preCount = input3d.extent(0);
	int postCount = input3d.extent(2);
	for (int i = 0; i < preCount; i++)
	{
		for (int j = 0; j < postCount; j++)
		{
			/*
			 * Copy wavefunction slice along bspline rank to temp
			 * array. We do this since LAPACK needs contiguous arrays.
			 */
			psiSlice = input3d(i, Range(fromStart, toEnd), j).copy();

			// Call on BSpline function to perform expansion
			output3d(i, Range(fromStart, toEnd) , j) = 
				bsplineObject->ExpandFunctionInBSplines(psiSlice);
		}
	}
}

/*
 *
 */
template<int Rank>
void BSplineTransform::InverseTransform(Wavefunction<Rank> &psi)
{
	BSpline::Ptr bsplineObject = 
		psi.GetRepresentation().GetRepresentation(BaseRank).GetBSplineObject();

	//Project input and output into 3D arrays with bsplineRank in the middle. 
	Array<cplx, 3> input3d = MapToRank3(psi.GetData(), BaseRank, 1);
	Array<cplx, 3> output3d = MapToRank3(psi.GetData(), BaseRank, 1);

	output3d = 0;
	int preCount = input.extent(0);
	int postCount = input.extent(2);
	for (int i=0; i<preCount; i++)
	{
		for (int j=0; j<postCount; j++)
		{
			output(i, thetaIndex, j) +=  legendre * input(i, lIndex, j);
		}
	}
}

} // Namespace

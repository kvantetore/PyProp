#ifndef BSPLINETRANSFORM_H
#define BSPLINETRANSFORM_H

#include "../../common.h"
#include "../../wavefunction.h"
#include "bspline.h"


namespace BSpline
{

/*
 * This class performs transformation between grid space
 * and "b-spline space". The grid is simply made up of the
 * quadrature points for the b-splines.
 */
template<int Rank>
class BSplineTransform
{
private:
	int BaseRank;
	int BSplineDataName;
	int BSplineGridDataName;
	BSpline::Ptr BSplineObject;
	blitz::Array<cplx, 1> TempData;

public:
	
	// Constructor
	BSplineTransform() :
		BaseRank(-1),
		BSplineDataName(-1),
		BSplineGridDataName(-1)
	{}
	~BSplineTransform() {}

	void SetupStep(Wavefunction<Rank> &psi, BSpline::Ptr bsplineObject, int baseRank);

	/*
	 * Transforms. We use functions from the BSPLINE class
	 * to expand wavefunction in b-splines or reconstruct
	 * it on the (quadrature) grid.
	 */
	void ForwardTransform(Wavefunction<Rank> &psi);
	void InverseTransform(Wavefunction<Rank> &psi);

	int GetBaseRank()
	{
		return BaseRank;
	}

	void SetBaseRank(int baseRank)
	{
		BaseRank = baseRank;
	}
};

} // Namespace

#endif


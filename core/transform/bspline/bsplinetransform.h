#ifndef BSPLINETRANSFORM_H
#define BSPLINETRANSFORM_H

#include "../../common.h"
#include "../../wavefunction.h"

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

public:
	
	// Constructor
	BSplineTransform() {}
	~BSplineTransform() {}

	void SetupStep(const Wavefunction<Rank> &psi, int baseRank);

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


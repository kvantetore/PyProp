#ifndef BSPLINEREPRESENTATION_H
#define BSPLINEREPRESENTATION_H

#include "../../common.h"
#include "../representation.h"
#include "../../utility/boostpythonhack.h"
#include "../../transform/bspline/bspline.h"
#include "../../utility/blitzblas.h"

namespace BSpline
{

/*
 * Represents the wavefunction in a b-spline basis
 */
class BSplineRepresentation : public Representation<1>
{

public:
	typedef blitz::Array<double, 1> VectorType;
	typedef shared_ptr<BSplineRepresentation> Ptr;

private:
	VectorType Grid;
	VectorType Weights;
	BSpline::Ptr BSplineObject;

	blitz::Array<double, 2> OverlapMatrixFullRow;
	blitz::Array<double, 2> OverlapMatrixFullCol;
public:
	
	//Constructors
	BSplineRepresentation() : OverlapMatrixFullRow(0), OverlapMatrixFullCol(0) {}
	virtual ~BSplineRepresentation() {}

	virtual Representation<1>::RepresentationPtr Copy()
	{
		return Representation<1>::RepresentationPtr(new BSplineRepresentation(*this));
	}

	/* ---------- Implementation of Representation<1> interface ----------- */
		
	//Returns the size of the grid 
	virtual blitz::TinyVector<int, 1> GetFullShape()
	{
		return Grid.extent(0);
	}

	/*
	 * Performs an inner product between two b-spline wavefunctions
	 * This should probably not be called to often, because a faster version
	 * will be available somewhere else
	 */
	virtual std::complex<double> InnerProduct(const Wavefunction<1>& w1, const Wavefunction<1>& w2)
	{
		throw std::runtime_error("Please do not call InnerProduct on the sub representations. Combined representation will take care of proper InnerProduct");
	}

	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetGlobalWeights(int rank)
	{
		throw std::runtime_error("BSplineRepresentation is a Non-Orthogonal representation. Use Overlap Matrix instead");
	}
	
	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong b-spline rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return Grid;
	}

	virtual OverlapMatrix::Ptr GetGlobalOverlapMatrix(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong b-spline rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}

		return BSplineObject->GetOverlap();		
	}

	virtual void ApplyConfigSection(const ConfigSection &config);
	void SetupRepresentation(BSpline::Ptr bsplineObject);

	/*
	 * Return b-spline object
	 */
	BSpline::Ptr GetBSplineObject() { return BSplineObject; }

};

typedef boost::shared_ptr<BSplineRepresentation> BSplineRepresentationPtr;

} //Namespace

#endif


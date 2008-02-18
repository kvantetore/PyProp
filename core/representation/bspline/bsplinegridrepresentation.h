#ifndef BSPLINEGRID_REPERESENTATION_H
#define BSPLINEGRID_REPRESENTATION_H

#include "../../common.h"
#include "../representation.h"
#include "../../utility/boostpythonhack.h"
#include "../../transform/bspline/bspline.h"

namespace BSpline
{

/*
 * Represents the wavefunction on b-spline quadrature 
 * grid points.
 */
class BSplineGridRepresentation : public Representation<1>
{
public:
	typedef blitz::Array<double, 1> VectorType;
	typedef shared_ptr<BSplineGridRepresentation> Ptr;

private:
	VectorType Grid;
	VectorType Weights;
	BSpline::Ptr BSplineObject;

public:
	//Constructors:
	BSplineGridRepresentation() {}
	virtual ~BSplineGridRepresentation() {}

	virtual Representation<1>::RepresentationPtr Copy()
	{
		return Representation<1>::RepresentationPtr(new BSplineGridRepresentation(*this));
	}

	//Returns the size of the grid
	virtual blitz::TinyVector<int, 1> GetFullShape()
	{
		return blitz::TinyVector<int, 1>();
	}

	/*
	 * Performs an inner product between two b-spline wavefunctions
	 * This should probably not be called to often, because a faster version
	 * will be available somewhere else
	 */
	virtual std::complex<double> InnerProduct(const Wavefunction<1>& w1, const Wavefunction<1>& w2)
	{
		throw std::runtime_error("BSplineGridRepresentation::InnerProduct is not implemented");
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

	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetLocalWeights(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong angular rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return this->GetDistributedModel()->GetLocalArray(Weights, rank);
	}

	void SetupRepresentation(BSpline::Ptr p);

	virtual void ApplyConfigSection(const ConfigSection &config);

};

typedef boost::shared_ptr<BSplineGridRepresentation> BSplineGridRepresentationPtr;

} //Namespace

#endif


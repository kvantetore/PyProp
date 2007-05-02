#ifndef REDUCEDSPHERICAL_THETAREPRESENTATION_H
#define REDUCEDSPHERICAL_THETAREPRESENTATION_H

#include "../../common.h"
#include "../representation.h"
#include "thetarange.h"

namespace ReducedSpherical
{

/** Represents the wavefunction in an angular (theta) basis
 * Where the wavefunction is symmetric with respect to phi
 * The distribution of (theta) can be chosen when creating the
 * omega range.
 */
typedef Representation<1> Representation1D;

class ThetaRepresentation : public Representation1D
{
public:
	ThetaRange Range;                          //The points chosen for the angular grid

	//Constructors:
	ThetaRepresentation() {}
	virtual ~ThetaRepresentation() {}

	void SetupRepresentation(int maxL);

	//Returns the size of the grid
	virtual blitz::TinyVector<int, 1> GetFullShape()
	{
		return blitz::TinyVector<int, 1>();
	}

	/*
	 * Performs an inner product between two radial wavefunctions
	 * This should probably not be called to often, because a faster version
	 * will be available in SphericalRepresentation
	 */
	virtual std::complex<double> InnerProduct(const Wavefunction<1>& w1, const Wavefunction<1>& w2)
	{
		throw std::runtime_error("ThetaRepresentation::InnerProduct is not implemented");
	}

	/** 
	Returns the portion of the grid local to the current processor.
	**/
	virtual blitz::Array<double, 1> GetGlobalGrid(int rank)
	{
		if (rank != GetBaseRank())
		{
			cout << "Warning: Trying to get the wrong angular rank. Got " << rank << ", expected " << GetBaseRank() <<  endl;
		}
		return Range.GetGrid();
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
		return this->GetDistributedModel().GetLocalArray(Range.GetWeights(), rank);
	}

	/** Apply config, and set up Range
	  */
	virtual void ApplyConfigSection(const ConfigSection &config);
};

typedef boost::shared_ptr<ThetaRepresentation> ThetaRepresentationPtr;

} //Namespace

#endif


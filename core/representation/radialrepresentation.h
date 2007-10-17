#ifndef RADIALREPRESENTATION_H
#define RADIALREPRESENTATION_H

#include "../common.h"
#include "cartesianrepresentation.h"

class RadialRepresentation : public CartesianRepresentation<1>
{
public:
	typedef shared_ptr<RadialRepresentation> Ptr;

	//Constructors
	RadialRepresentation() {}
	
	RadialRepresentation(CartesianRange &r0) :
		CartesianRepresentation<1>(r0)
	{ }

	virtual Representation<1>::RepresentationPtr Copy()
	{
		return Representation<1>::RepresentationPtr(new RadialRepresentation(*this));
	}

	virtual void ApplyConfigSection(const ConfigSection &cfg)
	{
		//Check that the rank specified in the config file is prop
		double rmax;
		int rcount;
		cfg.Get("rmax", rmax);
		cfg.Get("rcount", rcount);

		Range(0) = CartesianRange(-rmax, rmax, rcount);
	}
};

typedef boost::shared_ptr<RadialRepresentation> RadialRepresentationPtr;

#endif

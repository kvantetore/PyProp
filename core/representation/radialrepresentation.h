#ifndef RADIALREPRESENTATION_H
#define RADIALREPRESENTATION_H

#include "../common.h"
#include "cartesianrepresentation.h"

class RadialRepresentation : public CartesianRepresentation<1>
{
public:
	//Constructors
	RadialRepresentation() {}
	
	RadialRepresentation(CartesianRange &r0) :
		CartesianRepresentation<1>(r0)
	{ }
	
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

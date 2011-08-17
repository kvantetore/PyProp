#include "sphericalharmonicbasisrepresentation.h"

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

namespace SphericalBasis
{

void SphericalHarmonicBasisRepresentation::ApplyConfigSection(const ConfigSection &config)
{
	using namespace boost::python;
	
	object indexIterator = config.Get<object>("index_iterator");

	typedef stl_input_iterator<LmIndex> Iterator;
	Iterator begin(indexIterator);
	Iterator end;

	Range.BeginIndexList();
	for (Iterator i=begin; i!=end; i++)
	{
		Range.AddIndex(*i);
	}
	Range.EndIndexList();

	//cout << "Setup SphericalHarmonicBasisRepresentation complete (" << Range.Count() << ")" << endl;
}	

} //namespace


#include "coupledsphericalharmonicrepresentation.h"

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

namespace CoupledSpherical
{

void CoupledSphericalHarmonicRepresentation::ApplyConfigSection(const ConfigSection &config)
{
	using namespace boost::python;
	
	object indexIterator = config.Get<object>("index_iterator");

	typedef stl_input_iterator<CoupledIndex> Iterator;
	Iterator begin(indexIterator);
	Iterator end;

	Range.BeginIndexList();
	for (Iterator i=begin; i!=end; i++)
	{
		Range.AddIndex(*i);
	}
	Range.EndIndexList();

	cout << "Setup CoupledSphericalHarmonicRepresentation complete (" << Range.Count() << ")" << endl;
}	

} //namespace


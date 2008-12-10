#include "coupledsphericalharmonicrepresentation.h"

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

namespace CoupledSpherical
{

void CoupledSphericalHarmonicRepresentation::ApplyConfigSection(const ConfigSection &config)
{
	using namespace boost::python;
	
	object indexIterator = config.Get<object>("index_iterator");
	stl_input_iterator<tuple> begin(indexIterator);
	stl_input_iterator<tuple> end;

	Range.BeginIndexList();
	for (stl_input_iterator<tuple> i=begin; i!=end; i++)
	{
		tuple c = (*i);
		int l1 = extract<int>(c[0]);
		int l2 = extract<int>(c[1]);
		int L = extract<int>(c[2]);
		int M = extract<int>(c[3]);
		Range.AddIndex( CoupledIndex(l1, l2, L, M) );
	}
	Range.EndIndexList();

	cout << "Setup CoupledSphericalHarmonicRepresentation complete (" << Range.Count() << ")" << endl;
}	

} //namespace


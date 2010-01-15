#include <boost/python.hpp>
#include "redirect.h"

using namespace std;
using namespace boost::python;

class python_streambuf_impl 
{
private:
	object sys;
	
public:
	streambuf *oldbuf;

	python_streambuf_impl() : oldbuf(0)
	{
		sys = object(handle<>(PyImport_ImportModule("sys")));	
		object out = sys.attr("stdout");
	}

	object get_stdout()
	{
		return sys.attr("stdout");
	}
};

python_streambuf::python_streambuf() : pimpl_(new python_streambuf_impl) {}
python_streambuf::~python_streambuf() { delete pimpl_; pimpl_ = 0; } 

int python_streambuf::overflow(int c)
{
	char outc = traits_type::to_char_type(c);
	pimpl_->get_stdout().attr("write")(outc);
	return traits_type::not_eof(c);
}

std::streambuf* python_streambuf::get_oldbuf()
{
	return pimpl_->oldbuf;
}

void python_streambuf::set_oldbuf(std::streambuf* oldbuf)
{
	pimpl_->oldbuf = oldbuf;
}


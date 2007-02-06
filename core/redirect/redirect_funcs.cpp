#include "redirect.h"

python_streambuf* redirect_cout()
{
	python_streambuf *newbuf = new python_streambuf();
	std::streambuf* oldbuf = std::cout.rdbuf(newbuf);
	newbuf->set_oldbuf(oldbuf);

	return newbuf;
}

void restore_cout(python_streambuf *newbuf)
{
	std::streambuf* oldbuf = newbuf->get_oldbuf();
	std::cout.rdbuf(oldbuf);
}





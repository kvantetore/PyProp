#ifndef REDIRECT_H
#define REDIRECT_H

#include <core/common.h>
#include <streambuf>

class python_streambuf_impl;

class python_streambuf : public std::streambuf
{
private:
	 python_streambuf_impl* pimpl_;
	
public:
	python_streambuf();
	~python_streambuf();

	virtual int overflow(int c = EOF);

	std::streambuf* get_oldbuf();
	void set_oldbuf(std::streambuf* oldbuf);
};

#endif


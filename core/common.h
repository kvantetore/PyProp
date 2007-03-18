#ifndef COMMON_H
#define COMMON_H

#include <complex>
#include <cmath>
#include <blitz/array.h>
#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include "configuration.h"

#define sqr(x) ((x)*(x))

using boost::shared_ptr;
using std::cout;
using std::endl;

/*
Complex
*/
typedef std::complex<double> cplx;
const cplx I = cplx(0.0, 1.0);

/*
Equality operators for TinyVector
*/
template<class T, int Rank>
inline bool operator==(const blitz::TinyVector<T,Rank> &v1, const blitz::TinyVector<T,Rank> &v2)
{
	for (int i=0; i<Rank; i++)
	{
		if (v1(i) != v2(i))
		{
			return false;
		}
	}
	return true;
}

template<class T, int Rank>
inline bool operator!=(const blitz::TinyVector<T,Rank> &v1, const blitz::TinyVector<T,Rank> &v2)
{
	return !(v1 == v2);
}

/*
To/From string values
*/
template<class T> 
inline std::string ToString(const T &value)
{
	std::ostringstream strm;
	strm << value << std::flush;
	return strm.str();
}

template<class T>
inline T FromString(const std::string &value)
{
		throw std::runtime_error("Not implemented");
}

template<class T>
inline T max(const T& a, const T& b)
{
	return a > b ? a : b;
}

template<class T>
inline T min(const T& a, const T& b)
{
	return a < b ? a : b;
}

#ifdef PYPROP_ROUND
inline int round(double x)
{
	return static_cast<int>( (x > 0.0) ? (x + 0.5) : (x - 0.5) );
}
#endif

#endif


#ifndef LASERHELPER_H
#define LASERHELPER_H

class LaserHelper
{
public:
	static double C(double l, double m)
	{
		return l * std::sqrt( ((l+1+m)*(l+1-m)) / ((2*l+1)*(2*l+3)) );
	}

	static double D(double l, double m)
	{
		return -(l+1) * std::sqrt( ((l+m)*(l-m)) / ((2*l+1)*(2*l-1)) );
	}

	static double E(double l, double m)
	{
		return std::sqrt( ((l+1+m)*(l+1-m)) / ((2*l+1)*(2*l+3)) );
	}

	static double F(double l, double m)
	{
		return std::sqrt( ((l+m)*(l-m)) / ((2*l+1)*(2*l-1)) );
	}

	static int kronecker(int a, int b)
	{
		return a == b ? 1 : 0;
	}
};


#endif


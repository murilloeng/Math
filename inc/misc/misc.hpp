#pragma once

typedef void(*ndiff_fun)(double*, const double*, void**);

namespace math
{
	int sign(int);
	int sign(bool);
	int sign(double);

	void swap(double&, double&);
	void swap(unsigned&, unsigned&);

	double randu(double = 0, double = 1);
	double bound(double, double = -1, double = +1);

	void ndiff(ndiff_fun, double*, const double*, void**, unsigned, unsigned, double);
}
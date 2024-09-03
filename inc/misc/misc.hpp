#pragma once

//std
#include <cstdint>

typedef void(*ndiff_fun)(double*, const double*, void**);

namespace math
{
	int sign(int);
	int sign(bool);
	int sign(double);

	void swap(double&, double&);
	void swap(uint32_t&, uint32_t&);

	double randu(double = 0, double = 1);
	double bound(double, double = -1, double = +1);

	void ndiff(ndiff_fun, double*, const double*, void**, uint32_t, uint32_t, double);
}
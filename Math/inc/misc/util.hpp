#pragma once

//std
#include <ctime>
#include <cstdint>

typedef void(*ndiff_fun)(double*, const double*, void**);

namespace math
{
	int32_t sign(bool);
	int32_t sign(double);
	int32_t sign(int32_t);

	void swap(double&, double&);
	void swap(uint32_t&, uint32_t&);

	double randu(double = 0, double = 1);
	double bound(double, double = -1, double = +1);

	void fft(double*, const double*, uint32_t, bool);

	char* time_format(char*, const time_t&, bool = true);

	void ndiff(ndiff_fun, double*, const double*, void**, uint32_t, uint32_t, double);
}
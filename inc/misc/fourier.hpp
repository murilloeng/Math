#pragma once

//std
#include <cstdint>

namespace math
{
	void dft(double*, const double*, uint32_t);
	void idft(double*, const double*, uint32_t);
}
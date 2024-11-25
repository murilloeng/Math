//std
#include <cmath>

//math
#include "Math/inc/misc/fourier.hpp"

namespace math
{
	double fourier(const double* x, double t1, double t2, double w, uint32_t nt, uint32_t mode)
	{
		double z = 0;
		for(uint32_t i = 0; i < nt - 1; i++)
		{
			const double t = t1 + (t2 - t1) * (i + 0.5) / nt;
			z += (t2 - t1) / nt * (x[i] + x[i + 1]) / 2 * (mode ? sin(w * t) : cos(w * t));
		}
		return z;
	}
}
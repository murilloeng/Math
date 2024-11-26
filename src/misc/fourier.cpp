//std
#include <cmath>

//math
#include "Math/inc/misc/fourier.hpp"

namespace math
{
	void dft(double* z, const double* x, uint32_t n)
	{
		for(uint32_t i = 0; i < n; i++)
		{
			z[2 * i + 0] = 0;
			z[2 * i + 1] = 0;
			for(uint32_t j = 0; j < n; j++)
			{
				z[2 * i + 0] += x[j] * cos(2 * M_PI * i * j / n);
				z[2 * i + 1] -= x[j] * sin(2 * M_PI * i * j / n);
			}
		}
	}
	void idft(double* x, const double* z, uint32_t n)
	{
		for(uint32_t i = 0; i < n; i++)
		{
			x[i] = 0;
			for(uint32_t j = 0; j < n; j++)
			{
				x[i] += z[2 * j + 0] * cos(2 * M_PI * i * j / n);
				x[i] -= z[2 * j + 1] * sin(2 * M_PI * i * j / n);
			}
			x[i] /= n;
		}
	}
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
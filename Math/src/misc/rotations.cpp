//std
#include <cmath>

//math
#include "Math/Math/inc/misc/rotations.hpp"

namespace math
{
	double fn(double t, uint32_t n)
	{
		if(fabs(t) < 2 * M_PI)
		{
			return funt(t, n);
		}
		else
		{
			if(n % 2)
			{
				const uint32_t p = n / 2;
				return (p % 2 ? -1 : +1) / pow(t, n) * (sin(t) - sint(t, p));
			}
			else
			{
				const uint32_t p = n / 2;
				return (p % 2 ? -1 : +1) / pow(t, n) * (cos(t) - cost(t, p));
			}
		}
	}
	double dfn(double t, uint32_t n)
	{
		return t * (n * fn(t, n + 2) - fn(t, n + 1));
	}
	double funt(double t, uint32_t n)
	{
		//data
		int s = 1;
		uint32_t k = 0;
		double v = 0, z = 1, dv;
		double p = std::tgamma(n + 1);
		//compute
		while(true)
		{
			dv = s * z / p;
			if(v + dv == v)
			{
				return v;
			}
			else
			{
				k++;
				v += dv;
				s *= -1;
				z *= t * t;
				p *= 2 * k + n;
				p *= 2 * k + n - 1;
			}
		}
	}
	double cost(double t, uint32_t n)
	{
		int s = 1;
		uint32_t a = 1;
		double v = 0, p = 1;
		for(uint32_t k = 0; k < n; k++)
		{
			v += s * p / a;
			s *= -1;
			p *= t * t;
			a *= 2 * k + 1;
			a *= 2 * k + 2;
		}
		return v;
	}
	double sint(double t, uint32_t n)
	{
		int s = 1;
		uint32_t a = 1;
		double v = 0, p = t;
		for(uint32_t k = 0; k < n; k++)
		{
			v += s * p / a;
			s *= -1;
			p *= t * t;
			a *= 2 * k + 2;
			a *= 2 * k + 3;
		}
		return v;
	}
}
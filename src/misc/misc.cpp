//std
#include <cmath>
#include <cstring>
#include <algorithm>

//math
#include "Math/inc/misc/misc.hpp"

namespace math
{
	int sign(int x)
	{
		return x == 0 ? 0 : x < 0 ? -1 : +1;
	}
	int sign(bool t)
	{
		return t ? +1 : -1;
	}
	int sign(double x)
	{
		return x == 0 ? 0 : x < 0 ? -1 : +1;
	}

	void swap(double& a, double& b)
	{
		a = a + b;
		b = a - b;
		a = a - b;
	}
	void swap(unsigned& a, unsigned& b)
	{
		a ^= b;
		b ^= a;
		a ^= b;
	}

	double randu(double a, double b)
	{
		return a + rand() * (b - a) / RAND_MAX;
	}
	double bound(double v, double a, double b)
	{
		return fmax(fmin(v, b), a);
	}

	void ndiff(void(*fun)(double*, const double*), double* K, const double* x, unsigned nv, unsigned nx, double dx)
	{
		//data
		double* xp = (double*) alloca(nx * sizeof(double));
		double* f1 = (double*) alloca(nv * sizeof(double));
		double* f2 = (double*) alloca(nv * sizeof(double));
		double* f3 = (double*) alloca(nv * sizeof(double));
		double* f4 = (double*) alloca(nv * sizeof(double));
		//setup
		memcpy(xp, x, nx * sizeof(double));
		//derivative
		for(unsigned i = 0; i < nx; i++)
		{
			//1st state
			xp[i] -= dx;
			fun(f1, xp);
			//2nd state
			xp[i] -= dx;
			fun(f2, xp);
			//3rd state
			xp[i] += dx;
			xp[i] += dx;
			xp[i] += dx;
			fun(f3, xp);
			//4th state
			xp[i] += dx;
			fun(f4, xp);
			//derivative
			xp[i] -= dx;
			xp[i] -= dx;
			for(unsigned j = 0; j < nv; j++)
			{
				K[j + nv * i] = (8 * f3[j] - 8 * f1[j] + f2[j] - f4[j]) / 12 / dx;
			}
		}
	}
}
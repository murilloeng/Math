//std
#include <cmath>
#include <ctime>
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <malloc.h>
#include <stdexcept>

//ext
#include "external/cpp/inc/fftw3.h"

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/linear/matrix.hpp"

namespace math
{
	int32_t sign(int32_t x)
	{
		return x == 0 ? 0 : x < 0 ? -1 : +1;
	}
	int32_t sign(bool t)
	{
		return t ? +1 : -1;
	}
	int32_t sign(double x)
	{
		return x == 0 ? 0 : x < 0 ? -1 : +1;
	}

	void swap(double& a, double& b)
	{
		a = a + b;
		b = a - b;
		a = a - b;
	}
	void swap(uint32_t& a, uint32_t& b)
	{
		a ^= b;
		b ^= a;
		a ^= b;
	}

	void skip_lines(FILE* file, uint32_t lines)
	{
		char buffer[2048];
		for(uint32_t i = 0; i < lines; i++)
		{
			if(!fgets(buffer, sizeof(buffer), file))
			{
				throw std::runtime_error("Unable to skip line from file!");
			}
		}
	}

	double randu(double a, double b)
	{
		return a + rand() * (b - a) / RAND_MAX;
	}
	double bound(double v, double a, double b)
	{
		return fmax(fmin(v, b), a);
	}

	void fft(double* z, const double* x, uint32_t n, bool d)
	{
		//data
		fftw_complex* xd = (fftw_complex*) fftw_malloc(n * sizeof(fftw_complex));
		fftw_complex* zd = (fftw_complex*) fftw_malloc(n * sizeof(fftw_complex));
		fftw_plan plan = fftw_plan_dft_1d(n, xd, zd, d ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
		//setup
		for(uint32_t i = 0; i < n; i++)
		{
			xd[i][1] = 0;
			xd[i][0] = x[i];
		}
		//compute
		fftw_execute(plan);
		for(uint32_t i = 0; i < n; i++)
		{
			z[2 * i + 0] = zd[i][0];
			z[2 * i + 1] = zd[i][1];
		}
		//cleanup
		fftw_destroy_plan(plan);
		fftw_free(xd); fftw_free(zd);
	}

	char* time_format(char* string, const time_t& time, bool date)
	{
		if(!date)
		{
			const tm* c = localtime(&time);
			const char* format = "%02d:%02d:%02d";
			sprintf(string, format, c->tm_hour - 1, c->tm_min, c->tm_sec);
		}
		else
		{
			const tm* c = localtime(&time);
			const char* format = "%02d:%02d:%02d %02d/%02d/%04d";
			sprintf(string, format, c->tm_hour, c->tm_min, c->tm_sec, c->tm_mday, c->tm_mon + 1, c->tm_year + 1900);
		}
		return string;
	}

	void ndiff(ndiff_fun_1 fun, double* K, const double* x, void** a, uint32_t nv, uint32_t nx, double dx)
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
		for(uint32_t i = 0; i < nx; i++)
		{
			//1st state
			xp[i] -= dx;
			fun(f1, xp, a);
			//2nd state
			xp[i] -= dx;
			fun(f2, xp, a);
			//3rd state
			xp[i] += dx;
			xp[i] += dx;
			xp[i] += dx;
			fun(f3, xp, a);
			//4th state
			xp[i] += dx;
			fun(f4, xp, a);
			//derivative
			xp[i] -= dx;
			xp[i] -= dx;
			for(uint32_t j = 0; j < nv; j++)
			{
				K[j + nv * i] = (8 * f3[j] - 8 * f1[j] + f2[j] - f4[j]) / 12 / dx;
			}
		}
	}
	void ndiff(ndiff_fun_2 fun, double* K, const double* x, const void** a, uint32_t nv, uint32_t nx, double dx)
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
		for(uint32_t i = 0; i < nx; i++)
		{
			//1st state
			xp[i] -= dx;
			fun(f1, xp, a);
			//2nd state
			xp[i] -= dx;
			fun(f2, xp, a);
			//3rd state
			xp[i] += dx;
			xp[i] += dx;
			xp[i] += dx;
			fun(f3, xp, a);
			//4th state
			xp[i] += dx;
			fun(f4, xp, a);
			//derivative
			xp[i] -= dx;
			xp[i] -= dx;
			for(uint32_t j = 0; j < nv; j++)
			{
				K[j + nv * i] = (8 * f3[j] - 8 * f1[j] + f2[j] - f4[j]) / 12 / dx;
			}
		}
	}
	void ndiff(std::function<void(double*, const double*)> fun, double* K, const double* x, uint32_t nv, uint32_t nx, double dx)
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
		for(uint32_t i = 0; i < nx; i++)
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
			for(uint32_t j = 0; j < nv; j++)
			{
				K[j + nv * i] = (8 * f3[j] - 8 * f1[j] + f2[j] - f4[j]) / 12 / dx;
			}
		}
	}
}
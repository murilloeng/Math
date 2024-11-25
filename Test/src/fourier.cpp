//std
#include <cmath>
#include <vector>

//math
#include "Math/inc/misc/fourier.hpp"

//test
#include "Math/Test/inc/tests.hpp"

void test_fourier(void)
{
	//data
	const double t1 = 0;
	const double t2 = 80 * M_PI + 2 * M_PI / 10 * rand() / RAND_MAX;
	const double w1 = 0;
	const double w2 = 20;
	const uint32_t nt = 5000;
	const uint32_t nw = 5000;
	double* x = (double*) alloca(nt * sizeof(double));
	double* z = (double*) alloca(2 * nw * sizeof(double));
	//time
	for(uint32_t i = 0; i < nt; i++)
	{
		const double t = t1 + (t2 - t1) * i / nt;
		x[i] = sin(10 * t);
	}
	//fourier
	FILE* file = fopen("test.txt", "w");
	for(uint32_t i = 0; i < nw; i++)
	{
		const double w = w1 + (w2 - w1) * i / (nw - 1);
		z[2 * i + 0] = math::fourier(x, t1, t2, w, nt, 0);
		z[2 * i + 1] = math::fourier(x, t1, t2, w, nt, 1);
		fprintf(file, "%+.6e %+.6e %+.6e\n", w, z[2 * i + 0], z[2 * i + 1]);
	}
	fclose(file);
}
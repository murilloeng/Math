//std
#include <cmath>
#include <vector>
#include <cstdio>

//math
#include "Math/Math/inc/misc/util.hpp"

//test
#include "Math/Test/inc/miscellaneous.hpp"

void tests::miscellaneous::fft(void)
{
	//data
	const uint32_t n = 10000;
	const double T = 2 * M_PI * 100;
	FILE* file = fopen("test.txt", "w");
	double* x = (double*) alloca(n * sizeof(double));
	double* z = (double*) alloca(2 * n * sizeof(double));
	//fourier
	for(uint32_t i = 0; i < n; i++)
	{
		x[i] = 0;
		const double t = T * i / n;
		for(uint32_t j = 0; j < 10; j++)
		{
			x[i] += (j + 1) * cos(2 * (j + 1) * t);
		}
		
	}
	math::fft(z, x, n, true);
	//file
	for(uint32_t i = 0; i < n; i++)
	{
		fprintf(file, "%+.6e %+.6e %+.6e\n", 2 * M_PI * i / T, z[2 * i + 0], z[2 * i + 1]);
	}
	fclose(file);
}
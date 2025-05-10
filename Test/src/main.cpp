//std
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>

//math
#include "Math/Math/inc/linear/vector.hpp"

//test
#include "Math/Test/inc/tests.hpp"

int main(void)
{
	//test
	// tests::solvers::harmonic::pyramid();
	tests::solvers::newton_raphson::truss_von_mises();
	return EXIT_SUCCESS;

	//data
	bool test;
	math::matrix A(3, 3);
	math::vector r(3), dx(3);
	const double T = 8.00e+01;
	const double m = 1.00e+00;
	const double L = 1.00e+00;
	const double g = 9.81e+00;
	const double a = 5.00e-01;
	const double b = 2.50e-01;
	const double q = 0.40 * M_PI;
	const uint32_t step_max = 80000;
	const uint32_t iteration_max = 10;
	const double tolerance = 1.00e-05;
	FILE* file = fopen("data.txt", "w");

	//derived
	const double h = T / step_max;

	//data
	double l_old = 0, l_new = 0;
	double x1_old = +L * sin(q), v1_old = 0, a1_old;
	double x2_old = -L * cos(q), v2_old = 0, a2_old;
	double x1_new = +L * sin(q), v1_new = 0, a1_new;
	double x2_new = -L * cos(q), v2_new = 0, a2_new;

	//acceleration
	a1_new = a1_old = -2 * l_new * x1_new / m;
	a2_new = a2_old = -2 * l_new * x2_new / m - g;

	//time loop
	for(uint32_t step = 0; step < step_max; step++)
	{
		//predictor
		l_new = l_old;
		a1_new = a1_old;
		a2_new = a2_old;
		v1_new = v1_old + h * a1_old;
		v2_new = v2_old + h * a2_old;
		x1_new = x1_old + h * v1_old + h * h / 2 * a1_old;
		x2_new = x2_old + h * v2_old + h * h / 2 * a2_old;
		//iterations
		for(uint32_t iteration = 0; iteration < iteration_max; iteration++)
		{
			//residue
			r[0] = m * a1_new + 2 * x1_new * l_new;
			r[1] = m * a2_new + 2 * x2_new * l_new + m * g;
			r[2] = 2 * x1_new * a1_new + 2 * x2_new * a2_new;
			r[2] += 2 * v1_new * v1_new + 2 * v2_new * v2_new;
			test = r.norm() < tolerance;
			//check
			if(test)
			{
				l_old = l_new;
				x1_old = x1_new;
				x2_old = x2_new;
				v1_old = v1_new;
				v2_old = v2_new;
				a1_old = a1_new;
				a2_old = a2_new;
				printf("step: %3d iteration: %d l: %+.2e x1: %+.2e x2: %+.2e v1: %+.2e v2: %+.2e a1: %+.2e a2: %+.2e\n", step, iteration, l_new, x1_new, x2_new, v1_new, v2_new, a1_new, a2_new);
				fprintf(file, "%+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e\n", h * step, l_new, x1_new, x2_new, v1_new, v2_new, a1_new, a2_new);
				break;
			}
			//tangent
			A(0, 1) = A(1, 0) = A(2, 2) = 0;
			A(0, 2) = 2 * x1_new;
			A(1, 2) = 2 * x2_new;
			A(0, 0) = A(1, 1) = 2 * l_new + m / h / h / b;
			A(2, 0) = 2 * a1_new + 4 * v1_new * a / h / b + 2 * x1_new / h / h / b;
			A(2, 1) = 2 * a2_new + 4 * v2_new * a / h / b + 2 * x2_new / h / h / b;
			//update
			A.solve(dx, r);
			l_new -= dx[2];
			x1_new -= dx[0];
			x2_new -= dx[1];
			v1_new -= a / h / b * dx[0];
			v2_new -= a / h / b * dx[1];
			a1_new -= 1 / h / h / b * dx[0];
			a2_new -= 1 / h / h / b * dx[1];
		}
		if(!test)
		{
			printf("Solver failed!\n");
			break;
		}
	}
	fclose(file);
	//return
	return EXIT_SUCCESS;
}
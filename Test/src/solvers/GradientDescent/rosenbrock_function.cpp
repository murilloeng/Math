//std
#include <cstdio>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/GradientDescent.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::rosenbrock_function(void)
{
	//data
	math::solvers::GradientDescent solver;

	//solver
	solver.m_size = 2;
	solver.step_size(1.00e-3);
	solver.iteration_max(100000);
	solver.gradient([] (double* g, const double* x) {
		g[1] = 200 * (x[1] - x[0] * x[0]);
		g[0] = 2 * (x[0] - 1) + 400 * (x[0] * x[0] - x[1]) * x[0];
	});

	//setup
	solver.allocate();
	solver.m_x_old[0] = -1.20e+00;
	solver.m_x_old[1] = +1.20e+00;

	//solve
	solver.solve();

	//print
	printf("Status: %s\n", solver.m_status ? "✅" : "❌");
}
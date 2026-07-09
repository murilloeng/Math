//std
#include <cstdio>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/GradientDescent.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::double_quadratic(void)
{
	//data
	math::solvers::GradientDescent solver;

	//solver
	solver.size(2);
	solver.step_size(1.00e-02);
	solver.iteration_max(10000);
	solver.gradient([] (double* g, const double* x) {
		g[0] = 2 * x[0];
		g[1] = 4 * x[1];
	});

	//setup
	solver.allocate();
	solver.state_old(0, +1.00e+00);
	solver.state_old(1, +1.00e+00);

	//solve
	solver.solve();

	//print
	printf("Status: %s\n", solver.status() ? "✅" : "❌");
}
//std
#include <cstdio>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/GradientDescent.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::single_quartic(void)
{
	//data
	math::solvers::GradientDescent solver;

	//solver
	solver.size(1);
	solver.step_size(1.00e-02);
	solver.iteration_max(10000);
	solver.gradient([] (double* g, const double* x) { g[0] = 4 * x[0] * (x[0] * x[0] - 1); });

	//setup
	solver.allocate();
	solver.state_old(0, +1.00e-05);

	//solve
	solver.solve();

	//print
	printf("Status: %s\n", solver.status() ? "✅" : "❌");
}
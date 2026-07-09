//std
#include <cstdio>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/GradientDescent.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::himmelblau_function(void)
{
	//data
	math::solvers::GradientDescent solver;

	//solver
	solver.m_size = 2;
	solver.step_size(1.00e-02);
	solver.iteration_max(100);
	solver.gradient([] (double* g, const double* x) {
		g[0] = 4 * (x[0] * x[0] + x[1] - 11) * x[0] + 2 * (x[0] + x[1] * x[1] - 7);
		g[1] = 2 * (x[0] * x[0] + x[1] - 11) + 4 * (x[0] + x[1] * x[1] - 7) * x[1];
	});

	//setup
	solver.allocate();
	solver.m_x_old[0] = +3.00e+00;
	solver.m_x_old[1] = -1.00e+00;

	//solve
	solver.solve();

	//print
	printf("Status: %s\n", solver.m_status ? "✅" : "❌");
}
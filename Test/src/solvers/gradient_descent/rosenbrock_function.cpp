//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/gradient_descent.hpp"

//test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::rosenbrock_function(void)
{
	//data
	math::solvers::gradient_descent solver;

	//solver
	solver.m_size = 2;
	solver.m_silent = false;
	solver.m_step_size = 1e-3;
	solver.m_iteration_max = 100000;
	solver.m_gradient = [] (double* g, const double* x) {
		g[1] = 200 * (x[1] - x[0] * x[0]);
		g[0] = 2 * (x[0] - 1) + 400 * (x[0] * x[0] - x[1]) * x[0];
	};

	//setup
	solver.allocate();
	solver.m_x_new[0] = -1.2;
	solver.m_x_new[1] = +1.0;

	//solve
	solver.solve();

	//print
	math::vector(solver.m_x_new, solver.m_size).print("solution:");
}
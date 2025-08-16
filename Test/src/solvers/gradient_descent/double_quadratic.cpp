//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/gradient_descent.hpp"

//test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::double_quadratic(void)
{
	//data
	math::solvers::gradient_descent solver;

	//solver
	solver.m_size = 2;
	solver.m_silent = true;
	solver.m_step_size = 1e-2;
	solver.m_iteration_max = 10000;
	solver.m_gradient = [] (double* g, const double* x) {
		g[0] = 2 * x[0];
		g[1] = 4 * x[1];
	};

	//setup
	solver.allocate();
	solver.m_x_new[0] = 1;
	solver.m_x_new[1] = 1;

	//solve
	solver.solve();

	//print
	math::vector(solver.m_x_new, solver.m_size).print("solution:");
}
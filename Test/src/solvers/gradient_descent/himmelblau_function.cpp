//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/gradient_descent.hpp"

//test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::himmelblau_function(void)
{
	//data
	math::solvers::gradient_descent solver;

	//solver
	solver.m_size = 2;
	solver.m_silent = false;
	solver.m_step_size = 1e-2;
	solver.m_iteration_max = 10000;
	solver.m_gradient = [] (double* g, const double* x) {
		g[0] = 4 * (x[0] * x[0] + x[1] - 11) * x[0] + 2 * (x[0] + x[1] * x[1] - 7);
		g[1] = 2 * (x[0] * x[0] + x[1] - 11) + 4 * (x[0] + x[1] * x[1] - 7) * x[1];
	};

	//setup
	solver.allocate();
	solver.m_x_new[0] = 3;
	solver.m_x_new[1] = -1;

	//solve
	solver.solve();

	//print
	math::vector(solver.m_x_new, solver.m_size).print("solution:");
}
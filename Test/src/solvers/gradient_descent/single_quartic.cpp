//math
#include "Math/Math/inc/solvers/gradient_descent.hpp"

//test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::single_quartic(void)
{
	//data
	math::solvers::gradient_descent solver;

	//solver
	solver.m_size = 1;
	solver.m_step_size = 1e-2;
	solver.m_iteration_max = 10000;
	solver.m_gradient = [] (double* g, const double* x) { g[0] = 4 * x[0] * (x[0] * x[0] - 1); };

	//setup
	solver.allocate();
	solver.m_x_new[0] = -1e-5;

	//solve
	solver.solve();
}
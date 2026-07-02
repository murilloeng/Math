//Math
#include "Math/Test/inc/solvers.hpp"
#include "Math/inc/Solvers/GradientDescent.hpp"

void tests::solvers::gradient_descent::single_quartic(void)
{
	//data
	const double x0 = -1e-5;
	math::solvers::GradientDescent solver;

	//solver
	solver.size(1);
	solver.m_step_size = 1e-2;
	solver.iteration_max(10000);
	solver.m_gradient = [] (double* g, const double* x) { g[0] = 4 * x[0] * (x[0] * x[0] - 1); };

	//setup
	solver.allocate();
	solver.state_new(&x0);

	//solve
	solver.solve();
}
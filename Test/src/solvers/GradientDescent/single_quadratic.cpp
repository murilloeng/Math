//Math
#include "Math/inc/Solvers/GradientDescent.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::single_quadratic(void)
{
	//data
	const double x0 = 1;
	math::solvers::GradientDescent solver;

	//solver
	solver.size(1);
	solver.m_step_size = 1e-2;
	solver.iteration_max(10000);
	solver.m_gradient = [] (double* g, const double* x) { g[0] = x[0]; };

	//setup
	solver.allocate();
	solver.state_new(&x0);

	//solve
	solver.solve();
}
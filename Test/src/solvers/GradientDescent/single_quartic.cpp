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
	solver.m_size = 1;
	solver.m_step_size = 1e-2;
	solver.iteration_max(10000);
	solver.m_gradient = [] (double* g, const double* x) { g[0] = 4 * x[0] * (x[0] * x[0] - 1); };

	//setup
	solver.allocate();
	solver.m_x_old[0] = +1.00e-05;

	//solve
	solver.solve();

	//print
	double g;
	solver.m_gradient(&g, solver.m_x_new);
	math::Vector(&g, 1).print("Gradient:");
	math::Vector(solver.m_x_new, 1).print("Solution:");
}
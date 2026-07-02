//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/GradientDescent.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::himmelblau_function(void)
{
	//data
	const double x0[] = {3, -1};
	math::solvers::GradientDescent solver;

	//solver
	solver.m_silent = 2;
	solver.m_step_size = 1e-2;
	solver.m_iteration_max = 10000;
	solver.m_gradient = [] (double* g, const double* x) {
		g[0] = 4 * (x[0] * x[0] + x[1] - 11) * x[0] + 2 * (x[0] + x[1] * x[1] - 7);
		g[1] = 2 * (x[0] * x[0] + x[1] - 11) + 4 * (x[0] + x[1] * x[1] - 7) * x[1];
	};

	//setup
	solver.allocate();
	solver.m_x_old[0] = x0[0];
	solver.m_x_old[1] = x0[1];

	//solve
	solver.solve();

	//print
	math::Vector(solver.m_x_new, 2).print("solution:");
}
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
	solver.size(2);
	solver.silent(false);
	solver.m_step_size = 1e-2;
	solver.iteration_max(10000);
	solver.m_gradient = [] (double* g, const double* x) {
		g[0] = 4 * (x[0] * x[0] + x[1] - 11) * x[0] + 2 * (x[0] + x[1] * x[1] - 7);
		g[1] = 2 * (x[0] * x[0] + x[1] - 11) + 4 * (x[0] + x[1] * x[1] - 7) * x[1];
	};

	//setup
	solver.allocate();
	solver.state_new(x0);

	//solve
	solver.solve();

	//print
	math::Vector(solver.state_new(), 2).print("solution:");
}
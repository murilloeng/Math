//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/GradientDescent.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::rosenbrock_function(void)
{
	//data
	const double x0[] = {-1.2, +1.0};
	math::solvers::GradientDescent solver;

	//solver
	solver.size(2);
	solver.silent(false);
	solver.m_step_size = 1e-3;
	solver.iteration_max(100000);
	solver.m_gradient = [] (double* g, const double* x) {
		g[1] = 200 * (x[1] - x[0] * x[0]);
		g[0] = 2 * (x[0] - 1) + 400 * (x[0] * x[0] - x[1]) * x[0];
	};

	//setup
	solver.allocate();
	solver.state_new(x0);

	//solve
	solver.solve();

	//print
	math::Vector(solver.state_new(), 2).print("solution:");
}
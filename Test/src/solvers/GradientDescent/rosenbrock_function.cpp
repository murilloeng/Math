//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/GradientDescent.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::rosenbrock_function(void)
{
	//data
	math::solvers::GradientDescent solver;

	//solver
	solver.m_size = 2;
	solver.m_step_size = 1e-3;
	solver.iteration_max(100000);
	solver.m_gradient = [] (double* g, const double* x) {
		g[1] = 200 * (x[1] - x[0] * x[0]);
		g[0] = 2 * (x[0] - 1) + 400 * (x[0] * x[0] - x[1]) * x[0];
	};

	//setup
	solver.allocate();
	solver.m_x_old[0] = -1.20e+00;
	solver.m_x_old[1] = +1.20e+00;

	//solve
	solver.solve();

	//print
	double g[2];
	solver.m_gradient(g, solver.m_x_new);
	math::Vector(g, 2).print("Gradient:");
	math::Vector(solver.m_x_new, 2).print("Solution:");
}
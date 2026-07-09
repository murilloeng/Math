//std
#include <cstdio>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/GradientDescent.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::gradient_descent::single_quadratic(void)
{
	//data
	math::solvers::GradientDescent solver;
	
	//solver
	solver.m_size = 1;
	solver.step_size(1.00e-01);
	solver.iteration_max(1000);
	solver.gradient([] (double* g, const double* x) { g[0] = x[0]; });
	
	//setup
	solver.allocate();
	solver.m_x_old[0] = 1.00e+00;
	
	//solve
	solver.solve();
	
	//print
	printf("Status: %s\n", solver.m_status ? "✅" : "❌");
}
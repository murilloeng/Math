//std
#include <cmath>

//math
#include "Math/Math/inc/solvers/newmark.hpp"

//test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::newmark::single_pendulum(void)
{
	//data
	const double g = 9.81e+00;
	const double L = 1.00e+00;
	const double x0 = 0.00e+00;
	const double v0 = 0.9999 * 2 * sqrt(g / L);
	math::solvers::newmark solver;
	//setup
	solver.m_size = 1;
	solver.m_step_max = 2000;
	solver.m_t_max = 2.00e+01;
	solver.m_convergence.m_type = math::solvers::convergence::type::fixed;
	//initials
	solver.allocate();
	solver.m_x_new[0] = x0;
	solver.m_v_new[0] = v0;
	//forces
	solver.m_internal_force = [g, L](double* fi, const double* x, const double* v){
		fi[0] = g / L * sin(x[0]);
	};
	solver.m_external_force = [](double* fe, const double*, const double*, double t){
		fe[0] = 0;
	};
	//tangents
	solver.m_inertia = [](double* M, const double*){
		M[0] = 1;
	};
	solver.m_damping = [](double* C, const double*, const double*, double t){
		C[0] = 0;
	};
	solver.m_stiffness = [g, L](double* K, const double* x, const double*, const double*, double t){
		K[0] = g / L * cos(x[0]);
	};
	//solve
	solver.solve();
	//save
	solver.save("single-pendulum.txt");
}
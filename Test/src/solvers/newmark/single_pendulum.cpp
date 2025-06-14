#pragma once

//math
#include "Math/Math/inc/solvers/newmark.hpp"

//test
#include "Math/Test/inc/tests.hpp"

void tests::solvers::newmark::single_pendulum(void)
{
	//data
	const double x0 = 0.99 * M_PI;
	const double v0 = 0.00e+00;
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
	solver.m_internal_force = [](double* fi, const double* x, const double* v){
		fi[0] = sin(x[0]);
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
	solver.m_stiffness = [](double* K, const double* x, const double*, const double*, double t){
		K[0] = cos(x[0]);
	};
	//solve
	solver.solve();
	//save
	solver.save("single-pendulum.txt");
}
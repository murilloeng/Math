//std
#include <cmath>

//math
#include "Math/Math/inc/solvers/newmark.hpp"

//test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::newmark::single_dof(void)
{
	//data
	const uint32_t n = 40;
	const double m = 1.00e+00;
	const double c = 5.00e-02;
	const double k = 1.00e+00;
	const double f = 1.00e+00;
	const double w = 2.00e+00;
	const double x0 = 1.00e+00;
	const double v0 = 0.00e+00;
	const double w0 = sqrt(k / m);
	math::solvers::newmark solver;
	//setup
	solver.m_size = 1;
	solver.m_step_max = 2000;
	solver.m_t_max = 2 * M_PI * n / w0;
	solver.m_convergence.m_type = math::solvers::convergence::type::fixed;
	//initials
	solver.allocate();
	solver.m_x_new[0] = x0;
	solver.m_v_new[0] = v0;
	//forces
	solver.m_internal_force = [k, c](double* fi, const double* x, const double* v){
		fi[0] = k * x[0] + c * v[0];
	};
	solver.m_external_force = [f, w](double* fe, const double*, const double*, double t){
		fe[0] = f * sin(w * t);
	};
	//tangents
	solver.m_inertia = [m](double* M, const double*){
		M[0] = m;
	};
	solver.m_damping = [c](double* C, const double*, const double*, double t){
		C[0] = c;
	};
	solver.m_stiffness = [k](double* K, const double*, const double*, const double*, double t){
		K[0] = k;
	};
	//solve
	solver.solve();
	//save
	solver.save("single-dof.txt");
}
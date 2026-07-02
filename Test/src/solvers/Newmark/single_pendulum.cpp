//std
#include <cmath>

//Math
#include "Math/Test/inc/solvers.hpp"
#include "Math/inc/Solvers/Newmark.hpp"
#include "Math/inc/Validation/Validator.hpp"

void tests::solvers::newmark::single_pendulum(void)
{
	//data
	const double g = 9.81e+00;
	const double L = 1.00e+00;
	const double x0 = +M_PI_2;
	const double v0 = 0.00e+00;
	math::solvers::Newmark solver;
	math::validation::Validator validator;
	//setup
	solver.m_size = 1;
	solver.step_max(1000);
	solver.m_t_max = 1.00e+01;
	solver.convergence().type(math::solvers::Convergence::Type::Fixed);
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
	solver.save("Test/data/Solvers/Newmark/single pendulum/numeric.txt");
	//validator
	validator.create_item();
	validator.item(0)->load_numeric("Test/data/Solvers/Newmark/single pendulum/numeric.txt", 3, 0);
	validator.item(0)->load_reference("Test/data/Solvers/Newmark/single pendulum/reference.txt", 0, 1);
	//validation
	validator.validate();
}
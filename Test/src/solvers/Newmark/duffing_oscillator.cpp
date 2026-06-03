//std
#include <cmath>

//Math
#include "Math/Test/inc/solvers.hpp"
#include "Math/inc/Solvers/Newmark.hpp"
#include "Math/inc/Validation/Validator.hpp"

void tests::solvers::newmark::duffing_oscillator(void)
{
	//data
	const double np = 10;
	const double m = +1.00e+00;
	const double c = +5.00e-02;
	const double k = +1.00e+00;
	const double F = +5.00e-01;
	const double w = +2.00e+00;
	const double x0 = +0.00e+00;
	const double v0 = +0.00e+00;
	const double knl = 3.00e+00;
	math::solvers::Newmark solver;
	math::validation::Validator validator;
	//setup
	solver.m_size = 1;
	solver.m_step_max = 1000;
	solver.m_t_max = 2 * M_PI * np / sqrt(k / m);
	solver.m_convergence.m_type = math::solvers::Convergence::Type::Fixed;
	//initials
	solver.allocate();
	solver.m_x_new[0] = x0;
	solver.m_v_new[0] = v0;
	//forces
	solver.m_internal_force = [k, c, knl](double* fi, const double* x, const double* v){
		fi[0] = k * x[0] + c * v[0] + knl * x[0] * x[0] * x[0];
	};
	solver.m_external_force = [F, w](double* fe, const double*, const double*, double t){
		fe[0] = F * sin(w * t);
	};
	//tangents
	solver.m_inertia = [m](double* M, const double*){
		M[0] = m;
	};
	solver.m_damping = [c](double* C, const double*, const double*, double){
		C[0] = c;
	};
	solver.m_stiffness = [k, knl](double* K, const double* x, const double*, const double*, double){
		K[0] = k + 3 * knl * x[0] * x[0];
	};
	//solve
	solver.solve();
	//save
	solver.save("Test/data/Solvers/Newmark/duffing oscillator/numeric.txt");
	//validator
	validator.create_item();
	validator.item(0)->load_numeric("Test/data/Solvers/Newmark/duffing oscillator/numeric.txt", 3, 0);
	validator.item(0)->load_reference("Test/data/Solvers/Newmark/duffing oscillator/reference.txt", 0, 1);
	//validation
	validator.validate();
}
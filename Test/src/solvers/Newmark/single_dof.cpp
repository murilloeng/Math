//std
#include <cmath>

//Math
#include "Math/Test/inc/solvers.hpp"
#include "Math/inc/Solvers/Newmark.hpp"
#include "Math/inc/Validation/Validator.hpp"

void tests::solvers::newmark::single_dof(void)
{
	//data
	const uint32_t np = 10;
	const double m = 1.00e+00;
	const double c = 5.00e-02;
	const double k = 1.00e+00;
	const double f = 1.00e+00;
	const double w = 2.00e+00;
	const double x0 = 1.00e+00;
	const double v0 = 0.00e+00;
	const double w0 = sqrt(k / m);
	math::solvers::Newmark solver;
	math::validation::Validator validator;
	//setup
	solver.size(1);
	solver.step_max(1000);
	solver.time_max(2 * M_PI * np / w0);
	solver.convergence().type(math::solvers::Convergence::Type::Fixed);
	//initials
	solver.allocate();
	solver.state_old(0, x0);
	solver.velocity_old(0, v0);
	//forces
	solver.internal_force([k, c](double* fi, const double* x, const double* v){
		fi[0] = k * x[0] + c * v[0];
	});
	solver.external_force([f, w](double* fe, const double*, const double*, double t){
		fe[0] = f * sin(w * t);
	});
	//tangents
	solver.inertia([m](double* M, const double*){
		M[0] = m;
	});
	solver.damping([c](double* C, const double*, const double*, double t){
		C[0] = c;
	});
	solver.stiffness([k](double* K, const double*, const double*, const double*, double t){
		K[0] = k;
	});
	//solve
	solver.solve();
	//save
	solver.save("Test/data/Solvers/Newmark/single dof/numeric.txt");
	//validator
	validator.create_item();
	validator.item(0)->load_numeric("Test/data/Solvers/Newmark/single dof/numeric.txt", 3, 0);
	validator.item(0)->load_reference("Test/data/Solvers/Newmark/single dof/reference.txt", 0, 1);
	//validation
	validator.validate();
}
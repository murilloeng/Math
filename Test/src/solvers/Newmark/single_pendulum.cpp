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
	solver.size(1);
	solver.silent(true);
	solver.step_max(1000);
	solver.time_max(1.00e+01);
	solver.convergence().type_set(math::solvers::Convergence::Type::Fixed);
	//system
	solver.allocate();
	solver.state_old(0, x0);
	solver.velocity_old(0, v0);
	solver.system([g, L](double* fi, double* fe, double* K, double* C, double* M, const double* x, const double*, const double*, double){
		//tangent
		M[0] = 1;
		C[0] = 0;
		K[0] = g / L * cos(x[0]);
		//forces
		fe[0] = 0;
		fi[0] = g / L * sin(x[0]);
	});
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
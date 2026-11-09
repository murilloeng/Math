//std
#include <cmath>

//Math
#include "Math/Test/inc/solvers.hpp"
#include "Math/inc/Solvers/RungeKutta.hpp"
#include "Math/inc/Validation/Validator.hpp"

void tests::solvers::runge_kutta::duffing_oscillator(void)
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
	math::solvers::RungeKutta solver;
	math::validation::Validator validator;
	//setup
	solver.size(1);
	solver.silent(true);
	solver.step_max(1000);
	solver.time_max(2 * M_PI * np / sqrt(k / m));
	//system
	solver.allocate();
	solver.state_old(0, x0);
	solver.velocity_old(0, v0);
	solver.system([m, c, k, knl, F, w](double* M, double* fi, double* fe, const double* x, const double* v, double t){
		M[0] = m;
		fe[0] = F * sin(w * t);
		fi[0] = k * x[0] + c * v[0] + knl * x[0] * x[0] * x[0];
	});
	//solve
	solver.solve();
	//save
	solver.save("Test/data/Solvers/Runge Kutta/duffing oscillator/numeric.txt");
	//validator
	validator.create_item();
	validator.item(0)->load_numeric("Test/data/Solvers/Runge Kutta/duffing oscillator/numeric.txt", 3, 0);
	validator.item(0)->load_reference("Test/data/Solvers/Runge Kutta/duffing oscillator/reference.txt", 0, 1);
	//validation
	validator.validate();
}
//std
#include <cmath>

//Math
#include "Math/Test/inc/solvers.hpp"
#include "Math/inc/Solvers/RungeKutta.hpp"
#include "Math/inc/Validation/Validator.hpp"

void tests::solvers::runge_kutta::double_pendulum(void)
{
	//data
	const double g = 9.81e+00;
	const double m1 = 1.00e+00;
	const double m2 = 1.00e+00;
	const double l1 = 1.00e+00;
	const double l2 = 1.00e+00;
	math::solvers::RungeKutta solver;
	math::validation::Validator validator;
	//setup
	solver.size(2);
	solver.silent(true);
	solver.step_max(2000);
	solver.time_max(1.00e+01);
	//system
	solver.allocate();
	solver.velocity_old(0, 0);
	solver.velocity_old(1, 0);
	solver.state_old(0, +M_PI_4);
	solver.state_old(1, -M_PI_4);
	solver.system([g, m1, m2, l1, l2](double* M, double* fi, double* fe, const double* x, const double* v, double){
		//data
		const double q1 = x[0];
		const double q2 = x[1];
		const double v1 = v[0];
		const double v2 = v[1];
		const double dq = q2 - q1;
		const double mt = m1 + m2;
		//external force
		fe[0] = -mt * g * l1 * sin(q1);
		fe[1] = -m2 * g * l2 * sin(q2);
		//inertia
		M[0] = mt * l1 * l1;
		M[3] = m2 * l2 * l2;
		M[1] = M[2] = m2 * l1 * l2 * cos(dq);
		//internal force
		fi[0] = -m2 * l1 * l2 * sin(dq) * v2 * v2;
		fi[1] = +m2 * l1 * l2 * sin(dq) * v1 * v1;
	});
	//solve
	solver.solve();
	//save
	solver.save("Test/data/Solvers/Runge Kutta/double pendulum/numeric.txt");
	//validator
	validator.create_item();
	validator.create_item();
	validator.item(0)->load_numeric("Test/data/Solvers/Runge Kutta/double pendulum/numeric.txt", 6, 0);
	validator.item(1)->load_numeric("Test/data/Solvers/Runge Kutta/double pendulum/numeric.txt", 6, 3);
	validator.item(0)->load_reference("Test/data/Solvers/Runge Kutta/double pendulum/reference.txt", 0, 1);
	validator.item(1)->load_reference("Test/data/Solvers/Runge Kutta/double pendulum/reference.txt", 0, 2);
	//validation
	validator.validate();
}
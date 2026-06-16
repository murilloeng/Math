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
	const double q1 = +M_PI_4;
	const double q2 = -M_PI_4;
	const double v1 = +0.00e+00;
	const double v2 = +0.00e+00;
	math::solvers::RungeKutta solver;
	math::validation::Validator validator;
	//setup
	solver.m_size = 2;
	solver.m_step_max = 2000;
	solver.m_t_max = 1.00e+01;
	solver.m_convergence.m_type = math::solvers::Convergence::Type::Fixed;
	//initials
	solver.allocate();
	solver.m_x_new[0] = q1;
	solver.m_x_new[1] = q2;
	solver.m_v_new[0] = v1;
	solver.m_v_new[1] = v2;
	//forces
	solver.m_internal_force = [m2, l1, l2](double* fi, const double* x, const double* v){
		//data
		const double q1 = x[0];
		const double q2 = x[1];
		const double v1 = v[0];
		const double v2 = v[1];
		const double dq = q2 - q1;
		//force
		fi[0] = -m2 * l1 * l2 * sin(dq) * v2 * v2;
		fi[1] = +m2 * l1 * l2 * sin(dq) * v1 * v1;
	};
	solver.m_external_force = [g, m1, m2, l1, l2](double* fe, const double* x, const double*, double){
		//data
		const double q1 = x[0];
		const double q2 = x[1];
		const double mt = m1 + m2;
		//force
		fe[0] = -mt * g * l1 * sin(q1);
		fe[1] = -m2 * g * l2 * sin(q2);
	};
	//tangents
	solver.m_inertia = [m1, m2, l1, l2](double* M, const double* x){
		//data
		const double q1 = x[0];
		const double q2 = x[1];
		const double dq = q2 - q1;
		const double mt = m1 + m2;
		//inertia
		M[0] = mt * l1 * l1;
		M[3] = m2 * l2 * l2;
		M[1] = M[2] = m2 * l1 * l2 * cos(dq);
	};
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
//std
#include <cmath>
#include <vector>

//Math
#include "Math/inc/Solvers/Harmonic.hpp"
#include "Math/inc/Validation/Validator.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

static const uint32_t ns = 200;
static const double p = 5.00e-02;
static const double m = 1.00e+00;
static const double c = 5.00e-02;
static const double k = 1.00e+00;
static const double knl = 3.00e+00;
static const double dp0 = 1.00e-01;
static const double w_min = 6.00e-01;
static const double w_max = 1.50e+00;

static void inertia(double* M, const double*)
{
	M[0] = m;
}
static void damping(double* C, const double*, const double*)
{
	C[0] = c;
}
static void stiffness(double* K, const double* x, const double*, const double*, double, double, double)
{
	K[0] = k + 3 * knl * x[0] * x[0];
}

static void internal_force(double* fi, const double* x, const double* v)
{
	fi[0] = k * x[0] + knl * x[0] * x[0] * x[0] + c * v[0];
}
static void external_force(double* fe, const double*, double t, double w)
{
	fe[0] = cos(w * t);
}

//first harmonic balance, exact at each solution point:
//a1 = A^2 * (k + 3 / 4 * knl * A^2 - m * w^2) / p
//b1 = A^2 * c * w / p
static double reference(double x)
{
	return x;
}

void tests::solvers::harmonic::duffing(void)
{
	//data
	math::solvers::Harmonic solver;
	math::validation::Validator validator;
	//setup
	solver.load(p);
	solver.dofs(1);
	solver.silent(true);
	solver.step_max(ns);
	solver.harmonics(1);
	solver.watch_dof(1);
	solver.step_size(dp0);
	solver.frequency(w_min);
	solver.stop_criteria().load_max(w_max);
	solver.control(math::solvers::Harmonic::Control::Frequency);
	solver.continuation().type(math::solvers::Continuation::Type::ArcLengthCylindrical);
	solver.stop_criteria().add_type(math::solvers::StopCriteria::Type::LoadLimitMaximum);
	//system
	solver.inertia(inertia);
	solver.damping(damping);
	solver.stiffness(stiffness);
	solver.internal_force(internal_force);
	solver.external_force(external_force);
	//allocate
	solver.allocate();
	//solve
	solver.solve();
	//save
	solver.save("Test/data/Solvers/Harmonic/Duffing/numeric.txt");
	//validator
	const uint32_t rows = solver.step();
	std::vector<double> xa(rows), xb(rows), ya(rows), yb(rows);
	for(uint32_t i = 0; i < rows; i++)
	{
		ya[i] = solver.state_data(i, 1);
		yb[i] = solver.state_data(i, 2);
		const double w = solver.load_data(i);
		const double A2 = ya[i] * ya[i] + yb[i] * yb[i];
		xb[i] = A2 * c * w / p;
		xa[i] = A2 * (k + 0.75 * knl * A2 - m * w * w) / p;
	}
	validator.create_item();
	validator.create_item();
	validator.item(0)->tolerance(1.00e-05);
	validator.item(1)->tolerance(1.00e-05);
	validator.item(0)->function(reference);
	validator.item(1)->function(reference);
	validator.item(0)->load_numeric(xa.data(), ya.data(), rows);
	validator.item(1)->load_numeric(xb.data(), yb.data(), rows);
	//validation
	validator.validate();
}
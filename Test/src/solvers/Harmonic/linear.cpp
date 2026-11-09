//std
#include <cmath>

//Math
#include "Math/inc/Solvers/Harmonic.hpp"
#include "Math/inc/Validation/Validator.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

static const double k = 4.00e+00;
static const double c = 1.00e-01;
static const double m = 1.00e+00;
static const double p = 1.00e+00;

static double reference_cos(double w)
{
	const double a = k - m * w * w;
	return p * a / (a * a + c * c * w * w);
}
static double reference_sin(double w)
{
	const double a = k - m * w * w;
	return p * c * w / (a * a + c * c * w * w);
}

static void internal_force(double* fi, const double* x, const double* v)
{
	fi[0] = k * x[0] + c * v[0];
}
static void external_force(double* fe, const double*, double t, double w)
{
	fe[0] = cos(w * t);
}

static void inertia(double* M, const double*)
{
	M[0] = m;
}
static void damping(double* C, const double*, const double*)
{
	C[0] = c;
}
static void stiffness(double* K, const double*, const double*, const double*, double, double, double)
{
	K[0] = k;
}

void tests::solvers::harmonic::linear(void)
{
	//data
	const uint32_t ns = 1000;
	const double w0 = 1.00e-01;
	const double wf = 4.00e+00;
	math::solvers::Harmonic solver;
	math::validation::Validator validator;
	//setup
	solver.dofs(1);
	solver.load(p);
	solver.step_max(ns);
	solver.harmonics(1);
	solver.watch_dof(1);
	solver.frequency(1.00e-01);
	solver.step_size((wf - w0) / ns);
	solver.control(math::solvers::Harmonic::Control::Frequency);
	solver.continuation().type(math::solvers::Continuation::Type::LoadControl);
	//system
	solver.allocate();
	solver.inertia(inertia);
	solver.damping(damping);
	solver.stiffness(stiffness);
	solver.internal_force(internal_force);
	solver.external_force(external_force);
	//solve
	solver.solve();
	//save
	solver.save("Test/data/Solvers/Harmonic/Linear/numeric.txt");
	//validator
	validator.create_item();
	validator.create_item();
	validator.item(0)->tolerance(1.00e-05);
	validator.item(1)->tolerance(1.00e-05);
	validator.item(0)->function(reference_cos);
	validator.item(1)->function(reference_sin);
	validator.item(0)->load_numeric("Test/data/Solvers/Harmonic/Linear/numeric.txt", 3, 1);
	validator.item(1)->load_numeric("Test/data/Solvers/Harmonic/Linear/numeric.txt", 3, 2);
	//validation
	validator.validate();
}
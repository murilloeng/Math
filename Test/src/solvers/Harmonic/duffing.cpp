//std
#include <cmath>
#include <ctime>
#include <cstdio>

//Math
#include "Math/inc/Solvers/Harmonic.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

static const uint32_t ns = 200;
static const double l = 5.00e-02;
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

void tests::solvers::harmonic::duffing(void)
{
	//data
	math::solvers::Harmonic solver;
	//setup
	solver.load(l);
	solver.dofs(1);
	solver.step_max(ns);
	solver.harmonics(3);
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
}
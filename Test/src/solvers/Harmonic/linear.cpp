//std
#include <cmath>
#include <cstdio>

//Math
#include "Math/inc/Solvers/Harmonic.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

static const double k = 4.00e+00;
static const double c = 1.00e-01;
static const double m = 1.00e+00;

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
	//setup
	solver.dofs(1);
	solver.step_max(ns);
	solver.harmonics(1);
	solver.load(1.00e+00);
	solver.m_watch_dof = 1;
	solver.frequency(1.00e-01);
	solver.m_dp0 = (wf - w0) / ns;
	solver.control(math::solvers::Harmonic::Control::Frequency);
	solver.continuation().type(math::solvers::Continuation::Type::LoadControl);
	//system
	solver.inertia(inertia);
	solver.damping(damping);
	solver.stiffness(stiffness);
	solver.internal_force(internal_force);
	solver.external_force(external_force);
	//setup
	solver.allocate();
	//solve
	solver.solve();
	//save
	solver.save("Test/data/Solvers/Harmonic/Linear/numeric.txt");
}
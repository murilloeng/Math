//std
#include <cmath>
#include <cstdio>

//Math
#include "Math/inc/solvers/Harmonic.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

static const double k = 4.00e+00;
static const double c = 1.00e-01;
static const double m = 1.00e+00;

static void internal_force(double* fi, const double* x, const double* v)
{
	fi[0] = k * x[0] + c * v[0];
}
static void external_force(double* fe, const double*, const double*, double t, double w)
{
	fe[0] = cos(w * t);
}

static void inertia(double* M, const double*)
{
	M[0] = m;
}
static void damping(double* C, const double*, const double*, double)
{
	C[0] = c;
}
static void stiffness(double* K, const double*, const double*, const double*, double, double, double)
{
	K[0] = k;
}

void tests::solvers::harmonic::oscillator(void)
{
	//data
	math::solvers::Harmonic solver;
	//setup
	solver.m_dofs = 1;
	solver.m_l = 1.00e+00;
	solver.m_w = 1.00e-01;
	solver.m_harmonics = 1;
	solver.m_watch_dof = 1;
	solver.m_dp0 = 2.00e-10;
	solver.m_step_max = 1000;
	solver.m_quadrature_order = 20;
	solver.m_control = math::solvers::Harmonic::Control::Frequency;
	solver.m_continuation.m_type = math::solvers::Continuation::Type::ControlLoad;
	//system
	solver.m_inertia = inertia;
	solver.m_damping = damping;
	solver.m_stiffness = stiffness;
	solver.m_internal_force = internal_force;
	solver.m_external_force = external_force;
	//setup
	solver.allocate();
	//solve
	solver.solve();
	//save
	solver.save("oscillator.dat");
}
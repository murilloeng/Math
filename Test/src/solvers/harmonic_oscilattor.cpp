//std
#include <cmath>
#include <cstdio>

//math
#include "Math/Math/inc/solvers/harmonic.hpp"

//test
#include "Math/Test/inc/tests.hpp"

static const double k = 4.00e+00;
static const double c = 1.00e-01;
static const double m = 1.00e+00;

static void internal_force(double* fi, const double* d, const double* v, void** args)
{
	fi[0] = k * d[0] + c * v[0];
}
static void external_force(double* fe, double t, double w, const double* d, void** args)
{
	fe[0] = sin(w * t);
}
static void external_force_gradient(double* dfew, double t, double w, const double* d, void** args)
{
	dfew[0] = t * cos(w * t);
}

static void inertia(double* M, const double* d, void** args)
{
	M[0] = m;
}
static void damping(double* C, const double* d, const double* v, void** args)
{
	C[0] = c;
}
static void stiffness(double* K, double t, double w, double l, const double* d, const double* v, const double* a, void** args)
{
	K[0] = k;
}

void tests::solvers::harmonic_oscillator(void)
{
	//data
	math::harmonic solver;
	//setup
	solver.m_size = 1;
	solver.m_l_0 = 1.0;
	solver.m_w_0 = 0.1;
	solver.m_harmonics = 1;
	solver.m_dpg = 1.80e-02;
	solver.m_step_max = 100;
	solver.m_tolerance = 1e-5;
	solver.m_iteration_max = 10;
	solver.m_quadrature_order = 20;
	solver.m_parameter = math::harmonic_parameter::frequency;
	//system
	solver.m_inertia = inertia;
	solver.m_damping = damping;
	solver.m_stiffness = stiffness;
	solver.m_internal_force = internal_force;
	solver.m_external_force = external_force;
	solver.m_external_force_gradient = external_force_gradient;
	//solve
	solver.solve();
	//save
	solver.save();
}
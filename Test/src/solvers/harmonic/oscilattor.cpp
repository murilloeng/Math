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
	fe[0] = cos(w * t);
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

void tests::solvers::harmonic::oscillator(void)
{
	//data
	math::harmonic solver;
	//setup
	solver.m_size = 1;
	solver.m_l_0 = 1.0;
	solver.m_w_0 = 1.0;
	solver.m_harmonics = 1;
	solver.m_dpg = 2.00e-03;
	solver.m_step_max = 1000;
	solver.m_tolerance = 1e-5;
	solver.m_iteration_max = 10;
	solver.m_quadrature_order = 20;
	solver.m_control = math::harmonic_control::frequency;
	solver.m_strategy = math::harmonic_strategy::arc_length_spherical;
	//system
	solver.m_inertia = inertia;
	solver.m_damping = damping;
	solver.m_stiffness = stiffness;
	solver.m_internal_force = internal_force;
	solver.m_external_force = external_force;
	//solve
	if(solver.solve()) solver.save();
}
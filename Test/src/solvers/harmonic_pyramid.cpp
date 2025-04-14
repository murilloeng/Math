//std
#include <cmath>
#include <cstdio>

//math
#include "Math/Math/inc/solvers/harmonic.hpp"

//test
#include "Math/Test/inc/tests.hpp"

static const double ap = 8.00e+00;
static const double rp = 4.00e-02;
static const double bp = sqrt(1 + ap * ap);

static void inertia(double* M, const double* d, void** args)
{
	M[0] = 1;
}
static void damping(double* C, const double* d, const double* v, void** args)
{
	C[0] = 2 * rp;
}
static void stiffness(double* K, double t, double w, double l, const double* d, const double* v, const double* a, void** args)
{
	const double zp = 1 + d[0];
	const double lp = sqrt(zp * zp + ap * ap);
	const double sp = lp / bp;
	const double ep = log(sp);
	const double gep = 1 / sp;
	const double hep = -1 / sp / sp;
	K[0] = bp * bp * (gep * gep + ep * hep) * zp * zp / lp / lp;
	K[0] += bp * bp * ep * gep / sp * (1 - zp * zp / lp / lp);
}

static void internal_force(double* fi, const double* d, const double* v, void** args)
{
	const double zp = 1 + d[0];
	const double lp = sqrt(zp * zp + ap * ap);
	const double sp = lp / bp;
	const double ep = log(sp);
	const double gep = 1 / sp;
	fi[0] = bp * bp * bp * ep * gep * zp / lp + 2 * rp * v[0];
}
static void external_force(double* fe, double t, double w, const double* d, void** args)
{
	fe[0] = cos(w * t) / 2;
}

void tests::solvers::harmonic_pyramid(void)
{
	//data
	math::harmonic solver;
	//setup
	solver.m_size = 1;
	solver.m_l_0 = 0.06;
	solver.m_w_0 = 0.20;
	solver.m_harmonics = 2;
	solver.m_dpg = 1.60e-03;
	solver.m_step_max = 8000;
	solver.m_w_max = 1.80e+00;
	solver.m_iteration_max = 100;
	solver.m_tolerance = 1.00e-05;
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
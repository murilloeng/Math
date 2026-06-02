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
	solver.m_l = l;
	solver.m_dofs = 1;
	solver.m_w = w_min;
	solver.m_dp0 = dp0;
	solver.m_step_max = ns;
	solver.m_harmonics = 3;
	solver.m_watch_dof = 1;
	solver.m_stop_criteria.m_p_max = w_max;
	solver.m_control = math::solvers::Harmonic::Control::Frequency;
	solver.m_continuation.m_type = math::solvers::Continuation::Type::ArcLengthCylindrical;
	solver.m_stop_criteria.m_types |= uint32_t(math::solvers::StopCriteria::Type::LoadLimitMaximum);
	//system
	solver.m_inertia = inertia;
	solver.m_damping = damping;
	solver.m_stiffness = stiffness;
	solver.m_internal_force = internal_force;
	solver.m_external_force = external_force;
	//allocate
	solver.allocate();
	//solve
	solver.solve();
	//save
	solver.save("Test/data/solvers/harmonic/duffing/numeric.dat");
}
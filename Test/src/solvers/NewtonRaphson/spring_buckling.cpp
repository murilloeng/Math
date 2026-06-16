//std
#include <cmath>

//Math
#include "Math/Test/inc/solvers.hpp"
#include "Math/inc/Validation/Validator.hpp"
#include "Math/inc/Solvers/NewtonRaphson.hpp"

static const double Ks = 1.00e+00;
static const double Fs = 1.00e+00;
static const double Ls = 1.00e+00;
static const double qs = 1.00e-04;

static double function(double x)
{
	return Ks / (Fs * Ls) * x / sin(x + qs);
}

void tests::solvers::newton_raphson::spring_buckling(void)
{
	//data
	math::solvers::NewtonRaphson solver;
	math::validation::Validator validator;
	//setup
	solver.m_size = 1;
	solver.m_dp0 = 1.00e-02;
	solver.m_step_max = 400;
	solver.m_residue = [](double* r, double p, const double* x)
	{
		r[0] = p - Ks / (Fs * Ls) * x[0] / sin(x[0] + qs);
	};
	solver.m_tangent_1 = [](double* g, double p, const double* x)
	{
		g[0] = 1;
	};
	solver.m_tangent_2 = [](double* K, double p, const double* x)
	{
		K[0] = Ks / (Fs * Ls) * (1 - x[0] * cos(x[0] + qs) / sin(x[0] + qs)) / sin(x[0] + qs);
	};
	solver.m_continuation.m_type = math::solvers::Continuation::Type::ArcLengthSpherical;
	//setup
	solver.allocate();
	solver.m_p_new = 0;
	solver.m_x_new[0] = 0;
	//solve
	solver.solve();
	solver.save("Test/data/Solvers/Newton Raphson/Spring buckling/numeric.txt");
	//validator
	validator.create_item();
	validator.item(0)->function(function);
	validator.item(0)->load_numeric("Test/data/Solvers/Newton Raphson/Spring buckling/numeric.txt", 0, 1);
	//validation
	validator.validate();
}
//Math
#include "Math/Test/inc/solvers.hpp"
#include "Math/inc/Validation/Validator.hpp"
#include "Math/inc/Solvers/NewtonRaphson.hpp"

static double function(double x)
{
	return x * (1 - x * x);
}

void tests::solvers::newton_raphson::truss_von_mises(void)
{
	//data
	math::solvers::NewtonRaphson solver;
	math::validation::Validator validator;
	//setup
	solver.m_size = 1;
	solver.step_max(400);
	solver.m_dp0 = 1.00e-02;
	solver.m_residue = [](double* r, double p, const double* x)
	{
		r[0] = -p - x[0] * (x[0] * x[0] - 1);
	};
	solver.m_tangent_1 = [](double* g, double p, const double* x)
	{
		g[0] = -1;
	};
	solver.m_tangent_2 = [](double* K, double p, const double* x)
	{
		K[0] = 3 * x[0] * x[0] - 1;
	};
	solver.continuation().type(math::solvers::Continuation::Type::ArcLengthSpherical);
	//setup
	solver.allocate();
	solver.m_p_new = 0;
	solver.m_x_new[0] = 1;
	//solve
	solver.solve();
	solver.save("Test/data/Solvers/Newton Raphson/Truss von Mises/numeric.txt");
	//validator
	validator.create_item();
	validator.item(0)->function(function);
	validator.item(0)->load_numeric("Test/data/Solvers/Newton Raphson/Truss von Mises/numeric.txt", 0, 1);
	//validation
	validator.validate();
}
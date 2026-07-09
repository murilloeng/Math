//Math
#include "Math/Test/inc/solvers.hpp"
#include "Math/inc/Validation/Validator.hpp"
#include "Math/inc/Solvers/NewtonRaphson.hpp"

static double function(double x)
{
	return x * (1 - x * x);
}
static void residue(double* r, double p, const double* x)
{
	r[0] = -p - x[0] * (x[0] * x[0] - 1);
}
static void tangent_p(double* g, double p, const double* x)
{
	g[0] = -1;
}
static void tangent_x(double* K, double p, const double* x)
{
	K[0] = 3 * x[0] * x[0] - 1;
}

void tests::solvers::newton_raphson::truss_von_mises(void)
{
	//data
	math::solvers::NewtonRaphson solver;
	math::validation::Validator validator;
	//setup
	solver.size(1);
	solver.step_max(400);
	solver.step_size(1.00e-02);
	solver.continuation().type(math::solvers::Continuation::Type::ArcLengthSpherical);
	//system
	solver.residue(residue);
	solver.tangent_p(tangent_p);
	solver.tangent_x(tangent_x);
	//initial
	solver.allocate();
	solver.state_old(0, 1);
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
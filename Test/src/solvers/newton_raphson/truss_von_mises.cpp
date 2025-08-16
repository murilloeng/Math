//Math
#include "Math/Math/inc/solvers/newton_raphson.hpp"

//Test
#include "Math/Test/inc/solvers.hpp"

void tests::solvers::newton_raphson::truss_von_mises(void)
{
	//data
	math::solvers::newton_raphson solver;
	//setup
	solver.m_size = 1;
	solver.m_dp0 = 1.00e-02;
	solver.m_step_max = 400;
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
	solver.m_continuation.m_type = math::solvers::continuation::type::arc_length_spherical;
	//setup
	solver.allocate();
	solver.m_p_new = 0;
	solver.m_x_new[0] = 1;
	//solve
	solver.solve();
	solver.save("truss.txt");
}
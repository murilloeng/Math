//std
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include <arpack/arpack.h>

//Test
#include "Math/Test/inc/fem.hpp"
#include "Math/Test/inc/eigen.hpp"
#include "Math/Test/inc/groups.hpp"
#include "Math/Test/inc/solvers.hpp"
#include "Math/Test/inc/geometry.hpp"
#include "Math/Test/inc/rotations.hpp"
#include "Math/Test/inc/validation.hpp"
#include "Math/Test/inc/miscellaneous.hpp"

int main(void)
{
	try
	{
		// tests::eigen::dense_svd();
		// tests::eigen::dense_sym_std_full();
		// tests::eigen::dense_sym_gen_full();
		// tests::eigen::dense_sym_std_part();
		// tests::eigen::dense_sym_gen_part();
		// tests::eigen::dense_non_std_full();
		// tests::eigen::dense_non_std_full();

		// tests::eigen::sparse_svd();
		// tests::eigen::sparse_sym_std_part();
		// tests::eigen::sparse_sym_gen_part();
		// tests::eigen::sparse_non_std_part();
		tests::eigen::sparse_non_gen_part();

		// tests::solvers::newton_raphson::spring_buckling();
		// tests::solvers::newton_raphson::truss_von_mises();

		// tests::solvers::newmark::single_dof();
		// tests::solvers::newmark::single_pendulum();
		// tests::solvers::newmark::double_pendulum();
		// tests::solvers::newmark::duffing_oscillator();

		// tests::solvers::runge_kutta::single_dof();
		// tests::solvers::runge_kutta::single_pendulum();
		// tests::solvers::runge_kutta::double_pendulum();
		// tests::solvers::runge_kutta::duffing_oscillator();

		// tests::solvers::gradient_descent::single_quartic();
		// tests::solvers::gradient_descent::single_quadratic();
		// tests::solvers::gradient_descent::double_quadratic();
		// tests::solvers::gradient_descent::exponential_smooth();
		// tests::solvers::gradient_descent::himmelblau_function();
		// tests::solvers::gradient_descent::rosenbrock_function();
	}
	catch(const std::exception& exception)
	{
		printf("%s\n", exception.what());
	}
	//return
	return EXIT_SUCCESS;
}
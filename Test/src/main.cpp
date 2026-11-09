//std
#include <cstdio>
#include <stdexcept>

//Test
#include "Math/Test/inc/eigen.hpp"
#include "Math/Test/inc/groups.hpp"
#include "Math/Test/inc/solvers.hpp"
#include "Math/Test/inc/geometry.hpp"
#include "Math/Test/inc/validation.hpp"

int main(void)
{
	try
	{
		tests::solvers::harmonic::linear();
		tests::solvers::harmonic::duffing();

		tests::solvers::newmark::single_dof();
		tests::solvers::newmark::single_pendulum();
		tests::solvers::newmark::double_pendulum();
		tests::solvers::newmark::duffing_oscillator();

		tests::solvers::runge_kutta::single_dof();
		tests::solvers::runge_kutta::single_pendulum();
		tests::solvers::runge_kutta::double_pendulum();
		tests::solvers::runge_kutta::duffing_oscillator();

		tests::solvers::newton_raphson::spring_buckling();
		tests::solvers::newton_raphson::truss_von_mises();
	}
	catch(const std::exception& exception)
	{
		printf("%s\n", exception.what());
	}
	//return
	return EXIT_SUCCESS;
}
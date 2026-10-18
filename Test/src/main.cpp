//std
#include <cmath>
#include <cstdio>
#include <stdexcept>

//Test
#include "Math/inc/Groups/ASO2.hpp"
#include "Math/inc/Groups/GSO2.hpp"

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
		math::ASO2::test_adjoint(1000, 1e-5);
		math::ASO2::test_bracket(1000, 1e-5);
		math::ASO2::test_exp_log(1000, 1e-5);
		math::ASO2::test_jacobian(1000, 1e-5);
		math::ASO2::test_jacobian_inverse(1000, 1e-5);

		math::GSO2::test_matrix(1000, 1e-5);
		math::GSO2::test_inverse(1000, 1e-5);
		math::GSO2::test_log_exp(1000, 1e-5);
	}
	catch(const std::exception& exception)
	{
		printf("%s\n", exception.what());
	}
	//return
	return EXIT_SUCCESS;
}
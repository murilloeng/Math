//std
#include <cmath>
#include <cstdio>
#include <stdexcept>

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
		tests::eigen::non_symmetric_std();
		tests::eigen::non_symmetric_gen();
		tests::eigen::symmetric_std_full();
		tests::eigen::symmetric_gen_full();
		tests::eigen::symmetric_std_partial();
		tests::eigen::symmetric_gen_partial();
		tests::eigen::singular_value_decomposition();
	}
	catch(const std::exception& exception)
	{
		printf("%s\n", exception.what());
	}
	//return
	return EXIT_SUCCESS;
}
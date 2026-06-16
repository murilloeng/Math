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
		tests::solvers::harmonic::duffing();
	}
	catch(const std::exception& exception)
	{
		printf("%s\n", exception.what());
	}
	//return
	return EXIT_SUCCESS;
}
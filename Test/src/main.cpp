//std
#include <cstdio>
#include <stdexcept>

//test
#include "Math/Test/inc/fem.hpp"
#include "Math/Test/inc/groups.hpp"
#include "Math/Test/inc/solvers.hpp"
#include "Math/Test/inc/geometry.hpp"
#include "Math/Test/inc/rotations.hpp"
#include "Math/Test/inc/miscellaneous.hpp"

int main(void)
{
	try
	{
		//test
		tests::solvers::gradient_descent::exponential_smooth();
	}
	catch(const std::exception& exception)
	{
		printf("%s\n", exception.what());
	}
	//return
	return EXIT_SUCCESS;
}
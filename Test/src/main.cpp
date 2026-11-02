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
		tests::solvers::newton_raphson::spring_buckling();
	}
	catch(const std::exception& exception)
	{
		printf("%s\n", exception.what());
	}
	//return
	return EXIT_SUCCESS;
}
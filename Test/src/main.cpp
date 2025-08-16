//std
#include <ctime>
#include <cstdio>
#include <stdexcept>

//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/geometry/Circle.hpp"

//test
#include "Math/Test/inc/tests.hpp"

int main(void)
{
	try
	{
		//test
		tests::fem::beam::dynamics::rotation();
	}
	catch(const std::exception& exception)
	{
		printf("%s\n", exception.what());
	}
	//return
	return EXIT_SUCCESS;
}
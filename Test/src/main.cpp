//std
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdint>
#include <cstdlib>

//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/geometry/Circle.hpp"

//test
#include "Math/Test/inc/tests.hpp"

int main(void)
{
	//test
	tests::fem::beam::dynamics::section_strains_gradient();
	//return
	return EXIT_SUCCESS;
}
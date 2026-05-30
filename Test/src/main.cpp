//std
#include <cstdio>
#include <stdexcept>

//Test
#include "Math/Test/inc/fem.hpp"
#include "Math/Test/inc/groups.hpp"
#include "Math/Test/inc/solvers.hpp"
#include "Math/Test/inc/geometry.hpp"
#include "Math/Test/inc/rotations.hpp"
#include "Math/Test/inc/miscellaneous.hpp"

#include "Math/inc/solvers/Harmonic.hpp"

int main(void)
{
	try
	{
		math::solvers::Harmonic solver;
		// printf("state set: %d\n", solver.state_set());
		// tests::rotations::vec3::rotation_third();
	}
	catch(const std::exception& exception)
	{
		printf("%s\n", exception.what());
	}
	//return
	return EXIT_SUCCESS;
}
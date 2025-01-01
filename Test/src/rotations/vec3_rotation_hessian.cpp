//std
#include <ctime>
#include <cmath>
#include <cstdio>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat3.hpp"

//test
#include "Math/Test/inc/tests.hpp"

static bool inverse = true;
static bool transpose = true;

static void function(double* r, const double* v, void** args)
{
	//data
	math::vec3 rm = r;
	const math::vec3 vm = v;
	const math::vec3 am = (double*) args[0];
	//function
	rm = !inverse ? 
		vm.rotation_gradient(am, transpose) : 
		vm.rotation_gradient_inverse(am, transpose);
}
static void gradient(double* dr, const double* v, void** args)
{
	//data
	math::mat3 drm = dr;
	const math::vec3 vm = v;
	const math::vec3 am = (double*) args[0];
	//gradient
	drm = !inverse ? 
		vm.rotation_hessian(am, transpose) : 
		vm.rotation_hessian_inverse(am, transpose);
}

void tests::rotations::vec3_rotation_hessian(void)
{
	//data
	math::vec3 a, v, r;
	math::mat3 dra, drn, dri;
	const uint32_t nt = 10000;
	void* args[] = {a.data()};
	//test
	for(uint32_t i = 0; i < nt; i++)
	{
		a.randu();
		v.randu();
		bool test = true;
		gradient(dra.data(), v.data(), args);
		math::ndiff(function, drn.data(), v.data(), args, 3, 3, 1.00e-5);
		test = test && (dra - drn).norm() < 1e-5;
		printf("Test %04d: %s\n", i, test ? "ok" : "not ok");
		if(!test) break;
	}
}
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

static math::vec3 ar;

static bool inverse;
static bool transpose;

static void function(double* r, const double* v, void** args)
{
	//data
	math::vec3 rm = r;
	const math::vec3 vm = v;
	//function
	rm = !inverse ? 
		vm.rotation_gradient(ar, transpose) : 
		vm.rotation_gradient_inverse(ar, transpose);
}
static void gradient(double* dr, const double* v, void** args)
{
	//data
	math::mat3 drm = dr;
	const math::vec3 vm = v;
	//gradient
	drm = !inverse ? 
		vm.rotation_hessian(ar, transpose) : 
		vm.rotation_hessian_inverse(ar, transpose);
}

void tests::rotations::vec3_rotation_hessian(void)
{
	//data
	math::vec3 v, r;
	math::mat3 dra, drn, dri;
	const uint32_t nt = 10000;
	//menu
	while(true)
	{
		uint32_t selection;
		printf("Inverse?\n");
		printf("(1) Yes (2) No\n");
		scanf("%d", &selection);
		if(selection == 1 || selection == 2) {inverse = selection == 1; break;}
		printf("Invalid option!\n");
	}
	while(true)
	{
		uint32_t selection;
		printf("Transpose?\n");
		printf("(1) Yes (2) No\n");
		scanf("%d", &selection);
		if(selection == 1 || selection == 2) {transpose = selection == 1; break;}
		printf("Invalid option!\n");
	}
	//test
	for(uint32_t i = 0; i < nt; i++)
	{
		v.randu();
		ar.randu();
		bool test = true;
		gradient(dra.data(), v.data(), nullptr);
		math::ndiff(function, drn.data(), v.data(), nullptr, 3, 3, 1.00e-5);
		test = test && (dra - drn).norm() < 1e-5;
		printf("Test inverse(%d) transpose(%d) %04d: %s\n", inverse, transpose, i, test ? "ok" : "not ok");
		if(!test) break;
	}
}
//std
#include <ctime>
#include <cmath>
#include <cstdio>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/quat.hpp"
#include "Math/Math/inc/linear/mat3.hpp"

//test
#include "Math/Test/inc/rotations.hpp"

static bool mode;
static math::vec3 ar;
static math::quat qn;

/*
	r = exp(vs) * Rn * ar
	r = Rn * exp(vm) * ar
	dr = -spin(r) * T(vs) * dvs
	dr = -spin(r) * Rn * T(vm) * dvm
*/

static void function(double* r, const double* v, void** args)
{
	//data
	math::vec3 rm = r;
	const math::vec3 vm = v;
	//function
	rm = mode ? (vm.quaternion() * qn).rotate(ar) : (qn * vm.quaternion()).rotate(ar);
}
static void gradient(double* dr, const double* v, void** args)
{
	//data
	math::vec3 rm;
	math::mat3 drm = dr;
	const math::vec3 vm = v;
	function(rm.data(), v, args);
	//gradient
	if(mode)
	{
		drm = -rm.spin() * vm.rotation_gradient();
	}
	else
	{
		drm = -rm.spin() * qn.rotation() * vm.rotation_gradient();
	}
}

void tests::rotations::vec3::rotation_gradient(void)
{
	//data
	math::vec3 v, r;
	math::mat3 dra, drn, dri;
	const uint32_t nt = 10000;
	//menu
	while(true)
	{
		uint32_t selection;
		printf("Update mode:\n");
		printf("(1) Spatial (2) Material\n");
		const int args = scanf("%d", &selection);
		if(args == 1 && (selection == 1 || selection == 2)) {mode = selection == 1; break;}
		printf("Invalid option!\n");
	}
	//test
	for(uint32_t i = 0; i < nt; i++)
	{
		v.randu();
		ar.randu();
		qn.randu();
		bool test = true;
		dri = v.rotation_gradient_inverse();
		gradient(dra.data(), v.data(), nullptr);
		math::ndiff(function, drn.data(), v.data(), nullptr, 3, 3, 1.00e-5);
		test = test && (dra - drn).norm() < 1e-5;
		test = test && (v.rotation_gradient() * dri - math::mat3::eye()).norm() < 1e-5;
		printf("Test %s %d: %s\n", mode ? "spatial" : "material", i, test ? "ok" : "not ok");
		if(!test) break;
	}
}
//std
#include <ctime>
#include <cmath>
#include <cstdio>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat3.hpp"

//test
#include "Math/Test/inc/rotations.hpp"

static bool inverse;
static bool variable;
static bool transpose;
static math::vec3 t, u, v;

static void function(double* r, const double* x, void** args)
{
	//data
	math::vec3 rm = r;
	const math::vec3 tm = variable ? x : t;
	const math::vec3 am = variable ? u : x;
	//function
	rm = !inverse ? 
		tm.rotation_hessian(am, v, transpose) : 
		tm.rotation_hessian_inverse(am, v, transpose);
}
static void gradient(double* dr, const double* x, void** args)
{
	//data
	math::mat3 drm = dr;
	const math::vec3 tm = variable ? x : t;
	const math::vec3 am = variable ? u : x;
	//gradient
	drm = !inverse ? 
		tm.rotation_higher(
		am, v, transpose, variable) : 
		tm.rotation_higher(am, v, transpose, variable);
}

void tests::rotations::vec3::rotation_third(void)
{
	//data
	math::vec3 t, r;
	math::mat3 dra, drn, dri;
	const uint32_t nt = 10000;
	srand(uint32_t(time(nullptr)));
	const char* format = "Test inverse(%d) variable(%s) transpose(%d) %d: %s\n";
	//menu
	while(true)
	{
		uint32_t selection;
		printf("Inverse?\n");
		printf("(1) Yes (2) No\n");
		const int args = scanf("%d", &selection);
		if( args == 1 && (selection == 1 || selection == 2)) {inverse = selection == 1; break;}
		printf("Invalid option!\n");
	}
	while(true)
	{
		uint32_t selection;
		printf("Variable?\n");
		printf("(1) t (2) u\n");
		const int args = scanf("%d", &selection);
		if( args == 1 && (selection == 1 || selection == 2)) {variable = selection == 1; break;}
		printf("Invalid option!\n");
	}
	while(true)
	{
		uint32_t selection;
		printf("Transpose?\n");
		printf("(1) Yes (2) No\n");
		const int args = scanf("%d", &selection);
		if(args && (selection == 1 || selection == 2)) {transpose = selection == 1; break;}
		printf("Invalid option!\n");
	}
	//test
	for(uint32_t i = 0; i < nt; i++)
	{
		u.randu();
		v.randu();
		t.randu();
		bool test = true;
		gradient(dra.data(), t.data(), nullptr);
		math::ndiff(function, drn.data(), t.data(), nullptr, 3, 3, 1.00e-5);
		test = test && (dra - drn).norm() < 1e-5;
		printf(format, inverse, variable ? "t" : "u", transpose, i, test ? "ok" : "not ok");
		if(!test) break;
	}
}
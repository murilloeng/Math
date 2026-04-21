//std
#include <ctime>
#include <cmath>
#include <cstdio>

//math
#include "Math/inc/misc/util.hpp"
#include "Math/inc/linear/vec3.hpp"
#include "Math/inc/linear/mat3.hpp"

//test
#include "Math/Test/inc/rotations.hpp"

static bool coupled;
static bool inverse;
static bool transpose;
static math::vec3 u, v;

static void function(double* r, const double* t, void** args)
{
	//data
	const math::vec3 tm = t;
	const math::mat3 Tm = tm.rotation_gradient(transpose);
	const math::vec3 Tv = tm.rotation_gradient(v, transpose);
	const math::mat3 Tmi = tm.rotation_gradient_inverse(transpose);
	const math::vec3 Tvi = tm.rotation_gradient_inverse(v, transpose);
	//function
	r[0] = !inverse ? 
		coupled ? Tm.bilinear(u, v) : Tv.inner(u) : 
		coupled ? Tmi.bilinear(u, v) : Tvi.inner(u);
}
static void gradient(double* dr, const double* t, void** args)
{
	//data
	math::vec3 drm = dr;
	const math::vec3 tm = t;
	const math::mat3 H0 = tm.rotation_hessian(0, transpose);
	const math::mat3 H1 = tm.rotation_hessian(1, transpose);
	const math::mat3 H2 = tm.rotation_hessian(2, transpose);
	const math::mat3 Hv = tm.rotation_hessian(v, transpose);
	const math::mat3 Hi0 = tm.rotation_hessian_inverse(0, transpose);
	const math::mat3 Hi1 = tm.rotation_hessian_inverse(1, transpose);
	const math::mat3 Hi2 = tm.rotation_hessian_inverse(2, transpose);
	const math::mat3 Hiv = tm.rotation_hessian_inverse(v, transpose);
	//gradient
	if(!coupled)
	{
		drm = (!inverse ? Hv.transpose() : Hiv.transpose()) * u;
	}
	else
	{
		dr[0] = !inverse ? H0.bilinear(u, v) : Hi0.bilinear(u, v);
		dr[1] = !inverse ? H1.bilinear(u, v) : Hi1.bilinear(u, v);
		dr[2] = !inverse ? H2.bilinear(u, v) : Hi2.bilinear(u, v);
	}
}

void tests::rotations::vec3::rotation_hessian(void)
{
	//data
	math::vec3 t, r;
	uint32_t selection;
	math::mat3 dra, drn, dri;
	const uint32_t nt = 10000;
	//menu
	while(true)
	{
		printf("Mode?\n");
		printf("(1) Coupled (2) Uncoupled\n");
		const int args = scanf("%d", &selection);
		if(args == 1 && (selection == 1 || selection == 2)) {coupled = selection == 2; break;}
		printf("Invalid option!\n");
	}
	while(true)
	{
		printf("Inverse?\n");
		printf("(1) Yes (2) No\n");
		const int args = scanf("%d", &selection);
		if(args == 1 && (selection == 1 || selection == 2)) {inverse = selection == 1; break;}
		printf("Invalid option!\n");
	}
	while(true)
	{
		printf("Transpose?\n");
		printf("(1) Yes (2) No\n");
		const int args = scanf("%d", &selection);
		if(args && (selection == 1 || selection == 2)) {transpose = selection == 1; break;}
		printf("Invalid option!\n");
	}
	//test
	for(uint32_t i = 0; i < nt; i++)
	{
		t.randu();
		u.randu();
		v.randu();
		bool test = true;
		gradient(dra.data(), t.data(), nullptr);
		math::ndiff(function, drn.data(), t.data(), nullptr, 1, 3, 1.00e-5);
		test = test && (dra - drn).norm() < 1e-5;
		printf("Test - Mode: %s, Inverse: %s, Transpose: %s, Step: %04d, Status: %s\n", 
				coupled ? "Uncoupled" : "Coupled", inverse ? "Yes" : "No", transpose ? "Yes" : "No", i, test ? "ok" : "not ok");
		if(!test) break;
	}
}
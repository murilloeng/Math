//std
#include <ctime>
#include <cstdio>
#include <cstdlib>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat3.hpp"
#include "Math/Math/inc/linear/quat.hpp"

//test
#include "Math/Test/inc/fem.hpp"

static math::vec3 s1, ar;
static math::quat q1r, q2r;
static uint32_t nt = 10000;
static uint32_t what, order;
static math::vector d(4), ga(4), gn(4);
static math::matrix ha(4, 4), hn(4, 4);

/*
	uk = R1 * sk
	nk = R2 * sk
	R1 = exp(v1s) * R1n
	R1 = R1n * exp(v1m)
	R2 = exp(v2s) * R2n
	R2 = R2n * exp(v2m)
	R2 = exp(a * u1) * R1
	R2 = R1 * exp(a * s1)
	dt2s = dt1s + da * n1
	dt2m = R2^T * R1 * dt1m + da * s1
	dv2s = T^{-1}(v2s) * T(v1s) * dv1s + da * T^{-1}(v2s) * n1
*/

static void setup(void)
{
	d.randu();
	s1.randu();
	ar.randu();
	q1r.randu();
	q2r.randu();
}
static void state(double* v2sk, const double* d, void** args)
{
	const double a = d[3];
	const math::quat q1n = math::vec3(d).quaternion() * q1r;
	*v2sk = (q1n * (a * s1).quaternion() * q2r.conjugate()).pseudo().inner(ar);
}
static void gradient(double* dv2sk, const double* d, void** args)
{
	const double a = d[3];
	const math::vec3 v1s = d;
	const math::quat q1n = v1s.quaternion() * q1r;
	const math::quat q2n = q1n * (a * s1).quaternion();
	const math::vec3 v2s = (q2n * q2r.conjugate()).pseudo();
	dv2sk[3] = v2s.rotation_gradient_inverse(ar, true).inner(q2n.rotate(s1));
	math::vec3(dv2sk + 0) = v1s.rotation_gradient(true) * v2s.rotation_gradient_inverse(ar, true);
}
static void hessian(double* hv2sk, double* d, void** args)
{
	//data
	const double a = d[3];
	const math::vec3 v1s = d;
	const math::quat q1n = v1s.quaternion() * q1r;
	const math::quat q2n = q1n * (a * s1).quaternion();
	const math::vec3 v2s = (q2n * q2r.conjugate()).pseudo();
	//vector
	const math::vec3 n1 = q2n.rotate(s1);
	//gradient
	const math::mat3 T1 = v1s.rotation_gradient();
	const math::mat3 T2i = v2s.rotation_gradient_inverse();
	const math::vec3 t2i = v2s.rotation_gradient_inverse(ar, true);
	//hessian
	const math::mat3 H1 = v1s.rotation_hessian(t2i, true);
	const math::mat3 H2i = v2s.rotation_hessian_inverse(ar, true);
	//full hessian
	math::matrix(hv2sk, 4, 4).span(0, 3, 3, 1) = T1.transpose() * H2i * T2i * n1;
	math::matrix(hv2sk, 4, 4).span(0, 0, 3, 3) = T1.transpose() * H2i * T2i * T1 + H1;
	math::matrix(hv2sk, 4, 4).span(3, 3, 1, 1) = n1.transpose() * H2i * T2i * n1 - t2i.transpose() * n1.spin() * n1;
	math::matrix(hv2sk, 4, 4).span(3, 0, 1, 3) = n1.transpose() * H2i * T2i * T1 - t2i.transpose() * n1.spin() * T1;
}

static void menu_what(void)
{
	while(true)
	{
		printf("Test what:\n");
		printf("(1) Energy (2) Dependencies\n");
		const int args = scanf("%d", &what);
		if(args == 1 && (what == 1 || what == 2)) break;
		printf("Invalid option!\n");
	}
}
static void menu_order(void)
{
	while(true)
	{
		printf("Order:\n");
		printf("(1) Gradient (2) Hessian\n");
		const int args = scanf("%d", &order);
		if(args == 1 && (order == 1 || order == 2)) break;
		printf("Invalid option!\n");
	}
}

void test_gradient(void)
{
	srand((uint32_t) time(nullptr));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		gradient(ga.data(), d.data(), nullptr);
		math::ndiff(state, gn.data(), d.data(), nullptr, 1, 4, 1.00e-5);
		const bool test = (ga - gn).norm() < 1e-5;
		printf("Test %04d: %s\n", i, test ? "ok" : "not ok");
		if(!test)
		{
			ga.print("ga");
			gn.print("gn");
			(ga - gn).print("ge", 1e-5);
			break;
		}
	}
}
void test_hessian(void)
{
	srand((uint32_t) time(nullptr));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		hessian(ha.data(), d.data(), nullptr);
		math::ndiff(gradient, hn.data(), d.data(), nullptr, 4, 4, 1.00e-5);
		const bool test = (ha - hn).norm() < 1e-5;
		printf("Test %04d: %s\n", i, test ? "ok" : "not ok");
		if(!test)
		{
			ha.print("ha");
			hn.print("hn");
			(ha - hn).print("he", 1e-5);
			break;
		}
	}
}

void tests::fem::revolute_fixed(void)
{
	menu_what();
	menu_order();
	order == 1 ? test_gradient() : test_hessian();
}
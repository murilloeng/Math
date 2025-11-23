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

//data
static math::vector d(12);

static double lr;
static math::vec3 g0, w0;
static math::vec3 x10, x20;
static math::quat q10, q20;
static math::quat q1n, q2n;
static const uint32_t nt = 10000;

static void setup(void)
{
	//data
	d.randu();
	x10.randu();
	x20.randu();
	q10.randu();
	q20.randu();
	q1n.randu();
	q2n.randu();
	//length
	lr = q10.conjugate(q20).pseudo().rotation_gradient_inverse(q10.conjugate(x20 - x10)).norm();
	//strains
	w0 = q10.conjugate(q20).pseudo() / lr;
	g0 = (lr * w0).rotation_gradient_inverse(q10.conjugate(x20 - x10)) / lr - math::vec3(1, 0, 0);
}

static void strains(double* es, const double* d, void** args)
{
	//data
	math::vec3 g(es + 0), w(es + 3);
	const math::vec3 x1 = x10 + math::vec3(d + 0);
	const math::vec3 x2 = x20 + math::vec3(d + 6);
	const math::quat q1 = q10 * q1n * math::vec3(d + 3).quaternion();
	const math::quat q2 = q20 * q2n * math::vec3(d + 9).quaternion();
	//strains
	w = q1.conjugate(q2).pseudo() / lr;
	g = (lr * w).rotation_gradient_inverse(q1.conjugate(x2 - x1)) / lr - math::vec3(1, 0, 0);
}
static void strains_gradient(double* des, const double* d, void** args)
{
	//data
	math::matrix B(des, 6, 12, math::mode::zeros);
	const math::vec3 x1 = x10 + math::vec3(d + 0);
	const math::vec3 x2 = x20 + math::vec3(d + 6);
	const math::quat q1 = q10 * q1n * math::vec3(d + 3).quaternion();
	const math::quat q2 = q20 * q2n * math::vec3(d + 9).quaternion();
	//data
	const math::mat3 R1 = (q10 * q1n).rotation();
	const math::mat3 R2 = (q20 * q2n).rotation();
	const math::mat3 At = q1.conjugate().rotation();
	const math::mat3 T1 = math::vec3(d + 3).rotation_gradient();
	const math::mat3 T2 = math::vec3(d + 9).rotation_gradient();
	const math::mat3 Ti = q1.conjugate(q2).pseudo().rotation_gradient_inverse();
	const math::mat3 Hi = q1.conjugate(q2).pseudo().rotation_hessian_inverse(q1.conjugate(x2 - x1) / lr);
	//gradient
	B.span(0, 0, 3, 3) = -Ti * At / lr;
	B.span(0, 6, 3, 3) = +Ti * At / lr;
	B.span(3, 3, 3, 3) = -Ti * At * R1 * T1 / lr;
	B.span(3, 9, 3, 3) = +Ti * At * R2 * T2 / lr;
	B.span(0, 9, 3, 3) = +Hi * Ti * At * R2 * T2;
	B.span(0, 3, 3, 3) = -Hi * Ti * At * R1 * T1 + Ti * At * (x2 - x1).spin() * R1 * T1 / lr;
}
// static void strains_hessian(double* d2es, double* d, void** args)
// {
// }

// static void menu_what(void)
// {
// 	while(true)
// 	{
// 		printf("Test what:\n");
// 		printf("(1) Energy (2) Strains\n");
// 		const int args = scanf("%d", &what);
// 		if(args == 1 && (what == 1 || what == 2)) break;
// 		printf("Invalid option!\n");
// 	}
// }
// static void menu_order(void)
// {
// 	while(true)
// 	{
// 		printf("Order:\n");
// 		printf("(1) Gradient (2) Hessian\n");
// 		const int args = scanf("%d", &order);
// 		if(args == 1 && (order == 1 || order == 2)) break;
// 		printf("Invalid option!\n");
// 	}
// }

static void test_gradient(void)
{
	math::matrix Ba(6, 12);
	math::matrix Bn(6, 12);
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		strains_gradient(Ba.data(), d.data(), nullptr);
		math::ndiff(strains, Bn.data(), d.data(), nullptr, 6, 12, 1.00e-4);
		const bool test = (Ba - Bn).norm() < 1e-5;
		printf("Test %d: %s\n", i, test ? "ok" : "not ok");
		if(!test)
		{
			Ba.print("Ba");
			Bn.print("Bn");
			(Ba - Bn).print("Be", 1e-5);
			break;
		}
	}
}
// static void test_hessian(void)
// {
// 	srand((uint32_t) time(nullptr));
// 	for(uint32_t i = 0; i < nt; i++)
// 	{
// 		setup();
// 		hessian(ha.data(), d.data(), nullptr);
// 		math::ndiff(gradient, hn.data(), d.data(), nullptr, 4, 4, 1.00e-5);
// 		const bool test = (ha - hn).norm() < 1e-5;
// 		printf("Test %04d: %s\n", i, test ? "ok" : "not ok");
// 		if(!test)
// 		{
// 			ha.print("ha");
// 			hn.print("hn");
// 			(ha - hn).print("he", 1e-5);
// 			break;
// 		}
// 	}
// }

void tests::fem::beam3D(void)
{
	// menu_what();
	// menu_order();
	// order == 1 ? test_gradient() : test_hessian();
	test_gradient();
}
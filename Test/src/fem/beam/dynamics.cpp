//std
#include <ctime>
#include <cmath>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat3.hpp"
#include "Math/Math/inc/linear/quat.hpp"

//test
#include "Math/Test/inc/tests.hpp"

//data
static double l0;
static math::vec3 x01, x02;
static math::quat q1n, q2n;

static void setup_data(void)
{
	x01.randu();
	x02.randu();
	q1n.randu();
	q2n.randu();
	l0 = 0.1 + 0.9 * rand() / RAND_MAX;
}
static void strains(double* es, const double* d, void**)
{
	//data
	math::vec3 gs(es + 0);
	math::vec3 ws(es + 3);
	const math::vec3 u1(d + 0);
	const math::vec3 v1(d + 3);
	const math::vec3 u2(d + 6);
	const math::vec3 v2(d + 9);
	const math::vec3 e1(1, 0, 0);
	const math::vec3 x1 = x01 + u1;
	const math::vec3 x2 = x02 + u2;
	const math::quat q1 = q1n * v1.quaternion();
	const math::quat q2 = q2n * v2.quaternion();
	//strains
	ws = q1.conjugate(q2).pseudo() / l0;
	gs = (l0 * ws).rotation_gradient_inverse(q1.conjugate(x2 - x1)) / l0 - e1;
}
static void strains_gradient(double* Bs, const double* d, void**)
{
	//data
	double es[6];
	math::matrix Bm(Bs, 6, 12);
	const math::vec3 u1(d + 0);
	const math::vec3 v1(d + 3);
	const math::vec3 u2(d + 6);
	const math::vec3 v2(d + 9);
	const math::vec3 gs(es + 0);
	const math::vec3 ws(es + 3);
	const math::vec3 e1(1, 0, 0);
	const math::vec3 x1 = x01 + u1;
	const math::vec3 x2 = x02 + u2;
	const math::quat q1 = q1n * v1.quaternion();
	const math::quat q2 = q2n * v2.quaternion();
	//rotation
	strains(es, d, nullptr);
	const math::vec3 xr = q1.conjugate(x2 - x1) / l0;
	const math::mat3 T1 = v1.rotation_gradient(true);
	const math::mat3 T2 = v2.rotation_gradient(true);
	const math::mat3 Twi = (l0 * ws).rotation_gradient_inverse(true);
	const math::mat3 Awi = (l0 * ws).rotation_hessian_inverse(xr, false);
	//gradient
	Bm.zeros();
	Bm.span(3, 9) = +Twi * T2 / l0;
	Bm.span(3, 3) = -Twi.transpose() * T1 / l0;
	Bm.span(0, 0) = -Twi.transpose() * q1.conjugate().rotation() / l0;
	Bm.span(0, 6) = +Twi.transpose() * q1.conjugate().rotation() / l0;
	Bm.span(0, 9) = +Awi * Twi * T2;
	Bm.span(0, 3) = -Awi * Twi.transpose() * T1 + Twi.transpose() * xr.spin() * T1;
}

void tests::fem::beam::dynamics::section_strains_gradient(void)
{
	//data
	double d[12];
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	math::matrix Ba(6, 12), Bn(6, 12), Br(6, 12);
	//test
	setup_data();
	for(uint32_t i = 0; i < nt; i++)
	{
		//compute
		math::vector(d, 12).randu();
		strains_gradient(Ba.data(), d, nullptr);
		math::ndiff(strains, Bn.data(), d, nullptr, 6, 12, 1e-5);
		//check
		Br = Ba - Bn;
		if(Br.norm() > 1e-5)
		{
			Ba.print("Ba");
			Bn.print("Bn");
			Br.print("Br", 1e-5);
			break;
		}
		else
		{
			printf("Test strains gradient: %d ok!\n", i);
		}

	}
}
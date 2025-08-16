//std
#include <ctime>
#include <cmath>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat3.hpp"
#include "Math/Math/inc/linear/quat.hpp"

//test
#include "Math/Test/inc/fem.hpp"

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

static void compute_strains(double* es, const double* d)
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
	const math::quat q1 = v1.quaternion() * q1n;
	const math::quat q2 = v2.quaternion() * q2n;
	//strains
	ws = q1.conjugate(q2).pseudo() / l0;
	gs = (l0 * ws).rotation_gradient_inverse(q1.conjugate(x2 - x1)) / l0 - e1;
}
static void compute_strains_gradient(double* Bs, const double* d)
{
	//data
	double es[6];
	math::matrix Bm(Bs, 6, 12);
	const math::vec3 u1(d + 0);
	const math::vec3 v1(d + 3);
	const math::vec3 u2(d + 6);
	const math::vec3 v2(d + 9);
	const math::vec3 ws(es + 3);
	const math::vec3 x1 = x01 + u1;
	const math::vec3 x2 = x02 + u2;
	const math::quat q1 = v1.quaternion() * q1n;
	//rotation
	compute_strains(es, d);
	const math::mat3 T1 = v1.rotation_gradient();
	const math::mat3 T2 = v2.rotation_gradient();
	const math::vec3 xr = q1.conjugate(x2 - x1) / l0;
	const math::mat3 R1t = q1.conjugate().rotation();
	const math::mat3 Twi = (l0 * ws).rotation_gradient_inverse();
	const math::mat3 Awi = (l0 * ws).rotation_hessian_inverse(xr);
	//gradient
	Bm.zeros();
	Bm.span(0, 0) = -Twi * R1t / l0;
	Bm.span(0, 6) = +Twi * R1t / l0;
	Bm.span(3, 9) = +Twi * R1t * T2 / l0;
	Bm.span(3, 3) = -Twi * R1t * T1 / l0;
	Bm.span(0, 9) = +Awi * Twi * R1t * T2;
	Bm.span(0, 3) = -Awi * Twi * R1t * T1 + Twi * xr.spin() * R1t * T1;
}
void compute_rotation(double* R, double s, const double* d)
{
	//data
	double es[6];
	const math::vec3 v1(d + 3);
	const math::vec3 ws(es + 3);
	const math::quat q1 = v1.quaternion() * q1n;
	//rotation
	compute_strains(es, d);
	math::mat3(R + 0) = q1.rotation() * (s * ws).rotation_tensor();
}
void compute_rotation_gradient(double* B, double s, const double* d)
{
	//data
	double es[6];
	const math::vec3 v1(d + 3);
	const math::vec3 ws(es + 3);
	const math::mat3 I(math::mode::eye);
	const math::quat q1 = v1.quaternion() * q1n;
	const math::mat3 T1 = v1.rotation_gradient();
	const math::mat3 R1t = q1.conjugate().rotation();
	const math::mat3 Tws = (s * ws).rotation_gradient();
	const math::mat3 Twi = (l0 * ws).rotation_gradient_inverse();
	//gradient
	math::matrix(B, 3, 12).zeros();
	const math::mat3 Rwt = (-s * ws).rotation_tensor();
	math::matrix(B, 3, 12).span(0, 3) = Rwt * (I - s / l0 * Tws * Twi) * R1t * T1;
}

void test_strains(double* es, const double* d, const void**)
{
	compute_strains(es, d);
}
void test_rotation(double* r, const double* d, const void** args)
{
	//data
	math::mat3 R;
	const double s = *(double*) args[1];
	const math::vec3& a = *(math::vec3*) args[1];
	//rotation
	compute_rotation(R.data(), s, d);
	//return
	math::vec3(r + 0) = R * a;
}

void tests::fem::beam::dynamics::strains(void)
{
	//data
	double d[12];
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	math::matrix Ba(6, 12), Bn(6, 12), Br(6, 12);
	//test
	for(uint32_t i = 0; i < nt; i++)
	{
		//compute
		setup_data();
		math::vector(d, 12).randu();
		compute_strains_gradient(Ba.data(), d);
		math::ndiff(test_strains, Bn.data(), d, nullptr, 6, 12, 1e-5);
		//check
		Br = Ba - Bn;
		if(Br.norm() > 1e-3 * Ba.norm())
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
void tests::fem::beam::dynamics::rotation(void)
{
	//data
	math::mat3 R;
	math::vec3 a;
	double s, d[12];
	const uint32_t nt = 100000;
	const void* args[] = {&s, &a};
	srand(uint32_t(time(nullptr)));
	math::matrix Ba(3, 12), Bn(3, 12), Br(3, 12), Bd(3, 12);
	//test
	for(uint32_t i = 0; i < nt; i++)
	{
		//compute
		a.randu();
		setup_data();
		s = l0 * rand() / RAND_MAX;
		math::vector(d, 12).randu();
		compute_rotation_gradient(Bd.data(), s, d);
		math::ndiff(test_rotation, Bn.data(), d, args, 3, 12, 1e-9);
		// Ba = ((math::matrix&) (-R * a.spin())) * Bd;
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
			printf("Test rotation gradient: %d ok!\n", i);
		}
	}
}
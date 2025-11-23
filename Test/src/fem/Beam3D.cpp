//std
#include <ctime>
#include <cmath>
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
static math::vector es(6);
static math::vector ss(6);
static math::vector eh(9);
static math::vector sh(9);
static math::vector d(12);
static uint32_t what, order;
static math::matrix A(9, 6);
static math::matrix B(6, 12);
static math::matrix Kh(9, 9);

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
	Kh.randu();
	x10.randu();
	x20.randu();
	q10.randu();
	q20.randu();
	q1n.randu();
	q2n.randu();
	Kh = (Kh + Kh.transpose()) / 2;
	//length
	lr = q10.conjugate(q20).pseudo().rotation_gradient_inverse(q10.conjugate(x20 - x10)).norm();
	//strains
	w0 = q10.conjugate(q20).pseudo() / lr;
	g0 = (lr * w0).rotation_gradient_inverse(q10.conjugate(x20 - x10)) / lr - math::vec3(1, 0, 0);
}

static void local_es(void)
{
	//data
	math::vec3 g(es.data() + 0), w(es.data() + 3);
	const math::vec3 x1 = x10 + math::vec3(d.data() + 0);
	const math::vec3 x2 = x20 + math::vec3(d.data() + 6);
	const math::quat q1 = q10 * q1n * math::vec3(d.data() + 3).quaternion();
	const math::quat q2 = q20 * q2n * math::vec3(d.data() + 9).quaternion();
	//data
	const math::vec3 tr = q1.conjugate(q2).pseudo();
	const math::vec3 xr = q1.conjugate(x2 - x1) / lr;
	//strains
	w = tr / lr;
	g = tr.rotation_gradient_inverse(xr) - math::vec3(1, 0, 0);
}
static void local_eh(void)
{
	//data
	const math::vec3 e1(1, 0, 0);
	const math::vec3 e2(0, 1, 0);
	const math::vec3 e3(0, 0, 1);
	const math::vec3 g(es.data() + 0);
	const math::vec3 w(es.data() + 3);
	//strains
	eh[1] = g[1] - g0[1];
	eh[2] = g[2] - g0[2];
	eh[3] = w[0] - w0[0];
	eh[6] = w0[1] * w0[2] - w[1] * w[2];
	eh[0] = (e1 + g).inner(e1 + g) / 2 - (e1 + g0).inner(e1 + g0) / 2;
	eh[4] = (1 + g[0]) * w[1] - g[1] * w[0] - (1 + g0[0]) * w0[1] + g0[1] * w0[0];
	eh[5] = (1 + g[0]) * w[2] - g[2] * w[0] - (1 + g0[0]) * w0[2] + g0[2] * w0[0];
	eh[7] = (w[0] * w[0] + w[2] * w[2]) / 2 - (w0[0] * w0[0] + w0[2] * w0[2]) / 2;
	eh[8] = (w[0] * w[0] + w[1] * w[1]) / 2 - (w0[0] * w0[0] + w0[1] * w0[1]) / 2;
}
void local_sh(void)
{
	sh = Kh * eh;
}
void local_ss(void)
{
	//data
	const math::vec3 g(es.data() + 0);
	const math::vec3 w(es.data() + 3);
	//stresses
	ss[1] = sh[1] + g[1] * sh[0] - w[0] * sh[4];
	ss[2] = sh[2] + g[2] * sh[0] - w[0] * sh[5];
	ss[0] = (1 + g[0]) * sh[0] + w[1] * eh[4] + w[2] * eh[5];
	ss[4] = (1 + g[0]) * sh[4] + w[2] * sh[6] + w[1] * sh[8];
	ss[5] = (1 + g[0]) * sh[5] + w[1] * sh[6] + w[2] * sh[7];
	ss[3] = sh[3] - g[1] * sh[4] - g[2] * sh[5] + w[0] * (sh[7] + sh[8]);
}
void local_B(void)
{
	//data
	const math::vec3 x1 = x10 + math::vec3(d.data() + 0);
	const math::vec3 x2 = x20 + math::vec3(d.data() + 6);
	const math::quat q1 = q10 * q1n * math::vec3(d.data() + 3).quaternion();
	const math::quat q2 = q20 * q2n * math::vec3(d.data() + 9).quaternion();
	//data
	const math::vec3 tr = q1.conjugate(q2).pseudo();
	const math::vec3 xr = q1.conjugate(x2 - x1) / lr;
	//data
	const math::mat3 Xr = xr.spin();
	const math::mat3 R1 = (q10 * q1n).rotation();
	const math::mat3 R2 = (q20 * q2n).rotation();
	const math::mat3 At = q1.conjugate().rotation();
	const math::mat3 Ti = tr.rotation_gradient_inverse();
	const math::mat3 Hi = tr.rotation_hessian_inverse(xr);
	const math::mat3 T1 = math::vec3(d.data() + 3).rotation_gradient();
	const math::mat3 T2 = math::vec3(d.data() + 9).rotation_gradient();
	//gradient
	B.span(0, 0, 3, 3) = -Ti * At / lr;
	B.span(0, 6, 3, 3) = +Ti * At / lr;
	B.span(3, 3, 3, 3) = -Ti * At * R1 * T1 / lr;
	B.span(3, 9, 3, 3) = +Ti * At * R2 * T2 / lr;
	B.span(0, 9, 3, 3) = +Hi * Ti * At * R2 * T2;
	B.span(0, 3, 3, 3) = -Hi * Ti * At * R1 * T1 + Ti * Xr * At * R1 * T1;
}

static void strains(double* fun_es, const double* x, void** args)
{
	d = x;
	local_es();
	math::vector(fun_es, 6) = es;
}
static void strains_gradient(double* fun_des, const double* x, void** args)
{
	d = x;
	local_B();
	math::matrix(fun_des, 6, 12) = B;
}
static void strains_hessian(double* d2es, double* d, void** args)
{
	//data
	const math::vec3 f = ss.data() + 0;
	const math::vec3 m = ss.data() + 3;
	const math::vec3 x1 = x10 + math::vec3(d + 0);
	const math::vec3 x2 = x20 + math::vec3(d + 6);
	const math::quat q1 = q10 * q1n * math::vec3(d + 3).quaternion();
	const math::quat q2 = q20 * q2n * math::vec3(d + 9).quaternion();
	//data
	math::matrix H(d2es, 12, 12, math::mode::zeros);
	const math::vec3 tr = q1.conjugate(q2).pseudo();
	const math::vec3 xr = q1.conjugate(x2 - x1) / lr;
	//data
	const math::mat3 Xr = xr.spin();
	const math::mat3 A1 = q1.rotation();
	const math::mat3 R1 = (q10 * q1n).rotation();
	const math::mat3 R2 = (q20 * q2n).rotation();
	const math::mat3 At = q1.conjugate().rotation();
	const math::mat3 Ti = tr.rotation_gradient_inverse();
	const math::mat3 Hi = tr.rotation_hessian_inverse(xr);
	const math::mat3 Tt = tr.rotation_gradient_inverse(true);
	const math::mat3 T1 = math::vec3(d + 3).rotation_gradient();
	const math::mat3 T2 = math::vec3(d + 9).rotation_gradient();
	const math::vec3 g = m + lr * tr.rotation_hessian_inverse(xr).transpose() * f;
	//hessian
	// -q1.rotate(tr.rotation_gradient_inverse(f, true)); //0
	// +q1.rotate(tr.rotation_gradient_inverse(f, true)); //6
	// +q1.rotate(tr.rotation_gradient_inverse(g, true)); //9
	// -q1.rotate(tr.rotation_gradient_inverse(g, true)) - (x2 - x1).cross(f + 6); //3
	
	
	// -q1.rotate(Tt * f);

	// -A1 * Tt * f / lr; //0
	// (-Hi * Ti * At * R1 * T1 + Ti * Xr * At * R1 * T1).transpose() * f - (Ti * At * R1 * T1 / lr).transpose() * m; //3
	// +(Ti * At / lr).transpose() * f; //6
	// (+Hi * Ti * At * R2 * T2).transpose() * f + (Ti * At * R2 * T2 / lr).transpose() * m; //9
	
}
static void strains_hessian_function(double* r, const double* d, void** args)
{
	math::matrix B(6, 12);
	strains_gradient(B.data(), d, nullptr);
	math::vector(r, 12) = B.transpose() * ss;
}

static void energy(double* U, const double* x, void** args)
{
	d = x;
	local_es();
	local_eh();
	U[0] = Kh.bilinear(eh, eh) / 2;
}
static void energy_gradient(double* dU, const double* x, void** args)
{
	d = x;
	local_B();
	local_es();
	local_eh();
	local_sh();
	local_ss();
	math::vector(dU, 12) = B.transpose() * ss;
}
static void energy_hessian(double* d2U, const double* x, void** args)
{
	return;
}

static void menu_what(void)
{
	while(true)
	{
		printf("Test what:\n");
		printf("(1) Energy (2) Strains\n");
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

static void test_strains_gradient(void)
{
	math::matrix Ba(6, 12);
	math::matrix Bn(6, 12);
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		strains_gradient(Ba.data(), d.data(), nullptr);
		math::ndiff(strains, Bn.data(), d.data(), nullptr, 6, 12, 1.00e-5);
		const bool test = (Ba - Bn).norm() < 1e-5 * Ba.norm();
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
static void test_strains_hessian(void)
{
	math::vector r(12);
	math::matrix Ha(12, 12);
	math::matrix Hn(12, 12);
	math::matrix He(12, 12);
	srand((uint32_t) time(nullptr));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		strains_hessian(Ha.data(), d.data(), nullptr);
		math::ndiff(strains_hessian_function, Hn.data(), d.data(), nullptr, 12, 12, 1.00e-5);
		He = Ha - Hn;
		const bool test = He.norm() < 1e-5 * Ha.norm();
		printf("Test %04d: %s\n", i, test ? "ok" : "not ok");
		if(!test)
		{
			Ha.print("Ha");
			Hn.print("Hn");
			He.print("He", 1e-5);
			break;
		}
	}
}

static void test_energy_gradient(void)
{
	math::vector dUa(12);
	math::vector dUn(12);
	math::vector dUe(12);
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		energy_gradient(dUa.data(), d.data(), nullptr);
		math::ndiff(energy, dUn.data(), d.data(), nullptr, 1, 12, 1.00e-5);
		dUe = dUa - dUn;
		const bool test = dUe.norm() < 1e-5 * dUa.norm();
		printf("Test %d: %s\n", i, test ? "ok" : "not ok");
		if(!test)
		{
			dUa.print("dUa");
			dUn.print("dUn");
			dUe.print("dUe", 1e-5);
			break;
		}
	}
}
static void test_energy_hessian(void)
{
	math::matrix d2Ua(12, 12);
	math::matrix d2Un(12, 12);
	math::matrix d2Ue(12, 12);
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		energy_hessian(d2Ua.data(), d.data(), nullptr);
		math::ndiff(energy_gradient, d2Un.data(), d.data(), nullptr, 12, 12, 1.00e-5);
		d2Ue = d2Ua - d2Un;
		const bool test = d2Ue.norm() < 1e-5 * d2Ua.norm();
		printf("Test %d: %s\n", i, test ? "ok" : "not ok");
		if(!test)
		{
			d2Ua.print("d2Ua");
			d2Un.print("d2Un");
			d2Ue.print("d2Ue", 1e-5);
			break;
		}
	}
}

void tests::fem::beam3D(void)
{
	menu_what();
	menu_order();
	what == 1 ? 
		order == 1 ? test_energy_gradient() : test_energy_hessian() : 
		order == 1 ? test_strains_gradient() : test_strains_hessian();
}
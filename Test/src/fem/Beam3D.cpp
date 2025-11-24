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
static math::vector f(12);
static uint32_t what, order;
static math::matrix A(9, 6, math::mode::zeros);
static math::matrix B(6, 12, math::mode::zeros);
static math::matrix Ks(6, 6, math::mode::zeros);
static math::matrix Kh(9, 9, math::mode::zeros);
static math::matrix Kg(12, 12, math::mode::zeros);

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
	//data
	const math::vec3 xr0 = q10.conjugate(x20 - x10);
	const math::vec3 tr0 = q10.conjugate(q20).pseudo();
	//length
	lr = tr0.rotation_gradient_inverse(xr0).norm();
	//strains
	w0 = tr0 / lr;
	g0 = (lr * w0).rotation_gradient_inverse(xr0) / lr - math::vec3(1, 0, 0);
}

static void local_A(void)
{
	//data
	const math::vec3 g(es.data() + 0);
	const math::vec3 w(es.data() + 3);
	//matrix
	A(0, 1) = +g[1];
	A(0, 2) = +g[2];
	A(4, 3) = -g[1];
	A(5, 3) = -g[2];
	A(7, 3) = A(8, 3) = +w[0];
	A(4, 1) = A(5, 2) = -w[0];
	A(1, 1) = A(2, 2) = A(3, 3) = 1;
	A(4, 0) = A(6, 5) = A(8, 4) = w[1];
	A(5, 0) = A(6, 4) = A(7, 5) = w[2];
	A(0, 0) = A(4, 4) = A(5, 5) = 1 + g[0];
}
static void local_Ks(void)
{
	//material
	Ks = A.transpose() * Kh * A;
	//geometric
	Ks(0, 0) += sh[0];
	Ks(1, 1) += sh[0];
	Ks(2, 2) += sh[0];
	Ks(4, 4) += sh[8];
	Ks(5, 5) += sh[7];
	Ks(3, 3) += sh[7] + sh[8];
	Ks(0, 4) = Ks(4, 0) += sh[4];
	Ks(0, 5) = Ks(5, 0) += sh[5];
	Ks(4, 5) = Ks(5, 4) += sh[6];
	Ks(1, 3) = Ks(3, 1) -= sh[4];
	Ks(2, 3) = Ks(3, 2) -= sh[5];
}
static void local_B(void)
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
	const math::mat3 At = q1.conjugate().rotation();
	const math::mat3 Ti = tr.rotation_gradient_inverse();
	const math::mat3 Hi = tr.rotation_hessian_inverse(xr);
	//gradient
	B.span(0, 0, 3, 3) = -Ti * At / lr;
	B.span(0, 6, 3, 3) = +Ti * At / lr;
	B.span(3, 3, 3, 3) = -Ti * At / lr;
	B.span(3, 9, 3, 3) = +Ti * At / lr;
	B.span(0, 9, 3, 3) = +Hi * Ti * At;
	B.span(0, 3, 3, 3) = -Hi * Ti * At + Ti * Xr * At;
}
void local_dB(void)
{
	//data
	const math::vec3 fs(ss.data() + 0);
	const math::vec3 ms(ss.data() + 3);
	const math::vec3 x1 = x10 + math::vec3(d.data() + 0);
	const math::vec3 x2 = x20 + math::vec3(d.data() + 6);
	const math::quat q1 = q10 * q1n * math::vec3(d.data() + 3).quaternion();
	const math::quat q2 = q20 * q2n * math::vec3(d.data() + 9).quaternion();
	//data
	const math::vec3 tr = q1.conjugate(q2).pseudo();
	const math::vec3 xr = q1.conjugate(x2 - x1) / lr;
	//data
	const math::mat3 Xr = xr.spin();
	const math::mat3 A1 = q1.rotation();
	const math::mat3 Xd = (x2 - x1).spin();
	const math::mat3 At = q1.conjugate().rotation();
	const math::mat3 Ti = tr.rotation_gradient_inverse();
	const math::mat3 Tt = tr.rotation_gradient_inverse(true);
	const math::mat3 Hf = tr.rotation_hessian_inverse(fs, true);
	const math::mat3 Ht = tr.rotation_hessian_inverse(xr).transpose();
	const math::mat3 Pi = tr.rotation_higher_inverse(xr, fs, false, true);
	const math::mat3 Qi = tr.rotation_higher_inverse(xr, fs, false, false);
	const math::mat3 Hm = tr.rotation_hessian_inverse(ms / lr + Ht * fs, true);
	const math::mat3 F2 = +A1 * (Tt * fs).spin() * At / lr;
	//hessian
	Kg.span(6, 9) = +A1 * Hf * Ti * At / lr;
	Kg.span(9, 6) = +A1 * Tt * Qi * At / lr;
	Kg.span(9, 0) = -A1 * Tt * Qi * At / lr;
	Kg.span(9, 9) = +A1 * (Hm + Tt * Pi) * Ti * At;
	Kg.span(6, 3) = -A1 * (Hf * Ti + (Tt * fs).spin()) * At / lr;
	Kg.span(9, 3) = -A1 * ((Hm + Tt * Pi) * Ti - Tt * Qi * Xr + (Tt * (ms / lr + Ht * fs)).spin()) * At;
	//equilibrium
	Kg.span(0, 0) = -Kg.span3(6, 0);
	Kg.span(0, 3) = -Kg.span3(6, 3);
	Kg.span(0, 6) = -Kg.span3(6, 6);
	Kg.span(0, 9) = -Kg.span3(6, 9);
	Kg.span(3, 3) = -Kg.span3(9, 3) - Xd * Kg.span3(6, 3);
	Kg.span(3, 9) = -Kg.span3(9, 9) - Xd * Kg.span3(6, 9);
	Kg.span(3, 0) = -Kg.span3(9, 0) - Xd * Kg.span3(6, 0) - F2;
	Kg.span(3, 6) = -Kg.span3(9, 6) - Xd * Kg.span3(6, 6) + F2;
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
	eh[6] = w[1] * w[2] - w0[1] * w0[2];
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
	ss[0] = (1 + g[0]) * sh[0] + w[1] * sh[4] + w[2] * sh[5];
	ss[4] = (1 + g[0]) * sh[4] + w[2] * sh[6] + w[1] * sh[8];
	ss[5] = (1 + g[0]) * sh[5] + w[1] * sh[6] + w[2] * sh[7];
	ss[3] = sh[3] - g[1] * sh[4] - g[2] * sh[5] + w[0] * (sh[7] + sh[8]);
}

static void strains(double* fun_es, const double* x, void** args)
{
	d = x;
	local_es();
	math::vector(fun_es, 6) = es;
}
static void strains_gradient(double* fun_des, const double* x, void** args)
{
	//data
	d = x;
	const math::vec3 t1 = x + 3;
	const math::vec3 t2 = x + 9;
	math::matrix T(12, 12, math::mode::eye);
	const math::mat3 T1 = (q10 * q1n).rotation() * t1.rotation_gradient();
	const math::mat3 T2 = (q20 * q2n).rotation() * t2.rotation_gradient();
	//gradient
	local_B();
	T.span(3, 3) = T1;
	T.span(9, 9) = T2;
	math::matrix(fun_des, 6, 12) = B * T;
}
static void strains_hessian(double* d2es, double* x, void** args)
{
	//data
	d = x;
	local_B();
	const math::vec3 t1 = x + 3;
	const math::vec3 t2 = x + 9;
	math::matrix T(12, 12, math::mode::eye);
	const math::vector f = B.transpose() * ss;
	const math::mat3 T1 = (q10 * q1n).rotation() * t1.rotation_gradient();
	const math::mat3 T2 = (q20 * q2n).rotation() * t2.rotation_gradient();
	//hessian
	local_dB();
	T.span(3, 3) = T1;
	T.span(9, 9) = T2;
	const math::vec3 m1(f.data() + 3);
	const math::vec3 m2(f.data() + 9);
	math::matrix(d2es, 12, 12) = T.transpose() * Kg * T;
	math::matrix(d2es, 12, 12).span(3, 3) += t1.rotation_hessian((q10 * q1n).conjugate(m1), true);
	math::matrix(d2es, 12, 12).span(9, 9) += t2.rotation_hessian((q20 * q2n).conjugate(m2), true);
}
static void strains_hessian_function(double* r, const double* x, void** args)
{
	//data
	d = x;
	const math::vec3 t1 = x + 3;
	const math::vec3 t2 = x + 9;
	math::matrix T(12, 12, math::mode::eye);
	const math::mat3 T1 = (q10 * q1n).rotation() * t1.rotation_gradient();
	const math::mat3 T2 = (q20 * q2n).rotation() * t2.rotation_gradient();
	//gradient
	local_B();
	T.span(3, 3) = T1;
	T.span(9, 9) = T2;
	math::vector(r, 12) = T.transpose() * B.transpose() * ss;
}

static void energy(double* U, const double* x, void** args)
{
	d = x;
	local_es();
	local_eh();
	U[0] = lr * Kh.bilinear(eh, eh) / 2;
}
static void energy_gradient(double* dU, const double* x, void** args)
{
	//data
	d = x;
	const math::vec3 t1(x + 3);
	const math::vec3 t2(x + 9);
	math::matrix T(12, 12, math::mode::eye);
	T.span(3, 3) = (q10 * q1n).rotation() * t1.rotation_gradient();
	T.span(9, 9) = (q20 * q2n).rotation() * t2.rotation_gradient();
	//gradient
	local_B();
	local_es();
	local_eh();
	local_sh();
	local_ss();
	math::vector(dU, 12) = lr * T.transpose() * B.transpose() * ss;
}
static void energy_hessian(double* d2U, const double* x, void** args)
{
	//data
	d = x;
	const math::vec3 t1(x + 3);
	const math::vec3 t2(x + 9);
	math::matrix T(12, 12, math::mode::eye);
	T.span(3, 3) = (q10 * q1n).rotation() * t1.rotation_gradient();
	T.span(9, 9) = (q20 * q2n).rotation() * t2.rotation_gradient();
	//hessian
	local_B();
	local_es();
	local_A();
	local_eh();
	local_sh();
	local_ss();
	local_dB();
	local_Ks();
	const math::vector f = B.transpose() * ss;
	const math::vec3 m1(f.data() + 3);
	const math::vec3 m2(f.data() + 9);
	math::matrix(d2U, 12, 12) = lr * (T.transpose() * (B.transpose() * Ks * B + Kg) * T);
	math::matrix(d2U, 12, 12).span(3, 3) += lr * t1.rotation_hessian((q10 * q1n).conjugate(m1), true);
	math::matrix(d2U, 12, 12).span(9, 9) += lr * t2.rotation_hessian((q20 * q2n).conjugate(m2), true);
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
		ss.randu();
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

	// test_A();
	// return;
	menu_what();
	menu_order();
	what == 1 ? 
		order == 1 ? test_energy_gradient() : test_energy_hessian() : 
		order == 1 ? test_strains_gradient() : test_strains_hessian();
}
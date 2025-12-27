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
static double lr;
static uint32_t order;
static math::vec3 z1, z2;
static math::quat q1n, q2n;
static math::vec3 s1, s2, s3;
static const uint32_t nt = 10000;

static math::vector dl(6);
static math::vector fl(6);
static math::vector d(12);
static math::vector f(12);
static math::matrix Kl(6, 6);
static math::matrix T(6, 12, math::mode::zeros);

static void setup(void)
{
	//data
	d.randu();
	z1.randu();
	z2.randu();
	Kl.randu();
	q1n.randu();
	q2n.randu();
	Kl += Kl.transpose();
	//axes
	lr = (z2 - z1).norm();
	s1 = (z2 - z1) / lr;
	s1.triad(s2, s3);
}

static void compute_state_dl(void)
{
	//kinematics
	const math::vec3 x1 = z1 + math::vec3(d.data() + 0);
	const math::vec3 x2 = z2 + math::vec3(d.data() + 6);
	//nodal quaternions
	const math::quat q0(s1, s2, s3);
	const math::quat q1 = math::vec3(d.data() + 3).quaternion() * q1n;
	const math::quat q2 = math::vec3(d.data() + 9).quaternion() * q2n;
	//local state
	const math::quat qr = q1 * q0;
	dl.span(3, 0, 3, 1) = (qr.conjugate() * q2 * q0).pseudo();
	dl.span(0, 0, 3, 1) = qr.conjugate(x2 - x1) - math::vec3(lr, 0, 0);
}

static void compute_force_fl(void)
{
	fl = Kl * dl;
}

static void compute_kinematic_T(void)
{
	//data
	const math::vec3 u1 = d.data() + 0;
	const math::vec3 u2 = d.data() + 6;
	const math::vec3 t1 = d.data() + 3;
	const math::vec3 t2 = d.data() + 9;
	//positions
	const math::vec3 x1 = z1 + u1;
	const math::vec3 x2 = z2 + u2;
	//nodal quaternions
	const math::quat q0(s1, s2, s3);
	const math::quat q1 = t1.quaternion() * q1n;
	const math::quat q2 = t2.quaternion() * q2n;
	//gradient
	const math::quat qr = q1 * q0;
	const math::mat3 Xr = (x2 - x1).spin();
	const math::mat3 T1 = t1.rotation_gradient();
	const math::mat3 T2 = t2.rotation_gradient();
	const math::mat3 Rt = qr.conjugate().rotation();
	const math::vec3 tl = qr.conjugate(q2 * q0).pseudo();
	const math::mat3 Ti = tl.rotation_gradient_inverse();
	//gradient
	T.span(0, 0) = -Rt;
	T.span(0, 6) = +Rt;
	T.span(3, 3) = -Ti * Rt * T1;
	T.span(3, 9) = +Ti * Rt * T2;
	T.span(0, 3) = +Rt * Xr * T1;
}

static void internal_energy(double* U)
{
	compute_state_dl();
	*U = Kl.bilinear(dl) / 2;
}

static void internal_force(double* f)
{
	compute_state_dl();
	compute_force_fl();
	compute_kinematic_T();
	math::vector(f, 12) = T.transpose() * fl;
}

static void stiffness(double* K)
{
	//data
	const math::vec3 u1 = d.data() + 0;
	const math::vec3 u2 = d.data() + 6;
	const math::vec3 t1 = d.data() + 3;
	const math::vec3 t2 = d.data() + 9;
	//positions
	const math::vec3 x1 = z1 + u1;
	const math::vec3 x2 = z2 + u2;
	//nodal quaternions
	const math::quat q0(s1, s2, s3);
	const math::quat q1 = t1.quaternion() * q1n;
	const math::quat q2 = t2.quaternion() * q2n;
	//gradient
	const math::quat qr = q1 * q0;
	const math::mat3 Rr = qr.rotation();
	const math::mat3 Xr = (x2 - x1).spin();
	const math::mat3 T1 = t1.rotation_gradient();
	const math::mat3 T2 = t2.rotation_gradient();
	const math::mat3 Rt = qr.conjugate().rotation();
	const math::vec3 tl = qr.conjugate(q2 * q0).pseudo();
	const math::mat3 Ti = tl.rotation_gradient_inverse();
	//material part
	compute_state_dl();
	compute_force_fl();
	compute_kinematic_T();
	math::matrix(K, 12, 12) = T.transpose() * Kl * T;
	//geometric part
	const math::vec3 nl = fl.data() + 0;
	const math::vec3 ml = fl.data() + 3;
	const math::vec3 mp = Ti.transpose() * ml;
	const math::mat3 Hl = tl.rotation_hessian_inverse(ml, true);
	math::matrix(K, 12, 12).span(0, 3) += Rr * nl.spin() * Rt * T1;
	math::matrix(K, 12, 12).span(6, 3) -= Rr * nl.spin() * Rt * T1;
	math::matrix(K, 12, 12).span(9, 9) += t2.rotation_hessian(Rr * mp, true);
	math::matrix(K, 12, 12).span(3, 0) -= T1.transpose() * Rr * nl.spin() * Rt;
	math::matrix(K, 12, 12).span(3, 6) += T1.transpose() * Rr * nl.spin() * Rt;
	math::matrix(K, 12, 12).span(3, 3) += T1.transpose() * Rr * Hl * Ti * Rt * T1;
	math::matrix(K, 12, 12).span(9, 3) -= T2.transpose() * Rr * Hl * Ti * Rt * T1;
	math::matrix(K, 12, 12).span(3, 9) -= T1.transpose() * Rr * Hl * Ti * Rt * T2;
	math::matrix(K, 12, 12).span(9, 9) += T2.transpose() * Rr * Hl * Ti * Rt * T2;
	math::matrix(K, 12, 12).span(9, 3) -= T2.transpose() * Rr * mp.spin() * Rt * T1;
	math::matrix(K, 12, 12).span(3, 3) += T1.transpose() * Rr * mp.spin() * Rt * T1;
	math::matrix(K, 12, 12).span(3, 3) += T1.transpose() * Xr * Rr * nl.spin() * Rt * T1;
	math::matrix(K, 12, 12).span(3, 3) -= t1.rotation_hessian(Rr * mp + Xr * Rr * nl, true);
}

static void energy_function(double* U, const double* d, void** args)
{
	::d = d;
	internal_energy(U);
}
static void energy_gradient(double* f, const double* d, void** args)
{
	::d = d;
	internal_force(f);
}
static void energy_hessian(double* K, const double* d, void** args)
{
	::d = d;
	stiffness(K);
}

static void test_gradient(void)
{
	math::vector fa(12);
	math::vector fn(12);
	math::vector fe(12);
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		energy_gradient(fa.data(), d.data(), nullptr);
		math::ndiff(energy_function, fn.data(), d.data(), nullptr, 1, 12, 1.00e-5);
		fe = fa - fn;
		const bool test = fe.norm() < 1e-5 * fa.norm();
		printf("Test %d: %s\n", i, test ? "ok" : "not ok");
		if(!test)
		{
			fa.print("fa");
			fn.print("fn");
			fe.print("fe", 1e-5);
			break;
		}
	}
}
static void test_hessian(void)
{
	math::matrix Ka(12, 12);
	math::matrix Kn(12, 12);
	math::matrix Ke(12, 12);
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		energy_hessian(Ka.data(), d.data(), nullptr);
		math::ndiff(energy_gradient, Kn.data(), d.data(), nullptr, 12, 12, 1.00e-5);
		Ke = Ka - Kn;
		const bool test = Ke.norm() < 1e-5 * Ka.norm();
		printf("Test %d: %s\n", i, test ? "ok" : "not ok");
		if(!test)
		{
			Ka.print("Ka");
			Kn.print("Kn");
			Ke.print("Ke", 1e-5);
			break;
		}
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

void tests::fem::beam3DCR(void)
{
	menu_order();
	order == 1 ? test_gradient () : test_hessian();
}
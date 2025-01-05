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
#include "Math/Test/inc/tests.hpp"

static const uint32_t nt = 10000;
static const double s = 4.00e-02;
static const double ka = 1.00e+03;
static const double ks = 2.00e+03;
static const double kt = 3.00e+03;
static const double kb = 4.00e+03;

static math::quat q0, q1, q2;
static math::mat3 Ku, Kp, Kt, Tr, Nr;
static math::vec3 z1, z2, u1, u2, x1, x2, t1, t2, ur, xr, tr, fr, mp, mr, n1, n2, n3;

static void setup(void)
{
	//state
	z1.randu();
	z2.randu();
	n3.randu();
	n3 /= n3.norm();
	n3.triad(n1, n2);
	q0 = math::quat(n1, n2, n3);
	//stiffness
	Ku = {ks, 0, 0, 0, ks, 0, 0, 0, ka};
	Kp = {kb, 0, 0, 0, kb, 0, 0, 0, kt};
}
static void local_dof(void)
{
	//state
	xr = q1.conjugate(x2 - x1);
	ur = xr - s * math::vec3(1, 0, 0);
	tr = (q1.conjugate() * q2).pseudo();
	//force
	fr = Ku * ur;
	mp = Kp * tr;
	mr = tr.rotation_gradient_inverse(mp, true);
	//stiffness
	Tr = tr.rotation_gradient_inverse();
	Nr = tr.rotation_hessian_inverse(mp, true);
	Kt = Tr.transpose() * Kp * Tr + Nr * Tr;
}
static void global_dof(const double* d)
{
	//dof
	u1 = math::vec3(d + 0);
	t1 = math::vec3(d + 3);
	u2 = math::vec3(d + 6);
	t2 = math::vec3(d + 9);
	//state
	x1 = z1 + u1;
	x2 = z2 + u2;
	q1 = t1.quaternion() * q0;
	q2 = t2.quaternion() * q0;
}
static void energy(double* U, const double* d, void** args)
{
	global_dof(d);
	local_dof();
	U[0] = ur.inner(Ku * ur) / 2 + tr.inner(Kp * tr) / 2;
}
static void internal_force_1(double* fi, const double* d, void** args)
{
	global_dof(d);
	local_dof();
	math::vec3(fi + 0) = -q1.rotate(fr);
	math::vec3(fi + 6) = +q1.rotate(fr);
	math::vec3(fi + 9) = +q1.rotate(mr);
	math::vec3(fi + 3) = -q1.rotate(mr + xr.cross(fr));
}
static void internal_force_2(double* fi, const double* d, void** args)
{
	global_dof(d);
	local_dof();
	math::vec3(fi + 0) = -q1.rotate(fr);
	math::vec3(fi + 6) = +q1.rotate(fr);
	math::vec3(fi + 9) = +q1.rotate(mr);
	math::vec3(fi + 3) = -q1.rotate(mr + xr.cross(fr));
	math::vec3(fi + 3) = t1.rotation_gradient(math::vec3(fi + 3), true);
	math::vec3(fi + 9) = t2.rotation_gradient(math::vec3(fi + 9), true);
}
static void stiffness(double* K, const double* d, void** args)
{
	global_dof(d);
	local_dof();
	math::mat3 Sr = xr.spin();
	math::mat3 Fr = fr.spin();
	math::mat3 Mr = mr.spin();
	math::mat3 A1 = q1.rotation();
	for(uint32_t i = 0; i < 144; i++) K[i] = 1;
	math::matrix(K, 12, 12).span(0, 0) = +A1 * Ku * A1.transpose();
	math::matrix(K, 12, 12).span(6, 6) = +A1 * Ku * A1.transpose();
	math::matrix(K, 12, 12).span(0, 6) = -A1 * Ku * A1.transpose();
	math::matrix(K, 12, 12).span(6, 0) = -A1 * Ku * A1.transpose();
	math::matrix(K, 12, 12).span(9, 9) = +A1 * Kt * A1.transpose();
	math::matrix(K, 12, 12).span(3, 9) = -A1 * Kt * A1.transpose();
	math::matrix(K, 12, 12).span(0, 9) = math::mat3(math::mode::zeros);
	math::matrix(K, 12, 12).span(9, 0) = math::mat3(math::mode::zeros);
	math::matrix(K, 12, 12).span(6, 9) = math::mat3(math::mode::zeros);
	math::matrix(K, 12, 12).span(9, 6) = math::mat3(math::mode::zeros);
	math::matrix(K, 12, 12).span(9, 3) = -A1 * (Kt + Mr) * A1.transpose();
	math::matrix(K, 12, 12).span(3, 0) = +A1 * (Sr * Ku - Fr) * A1.transpose();
	math::matrix(K, 12, 12).span(3, 6) = -A1 * (Sr * Ku - Fr) * A1.transpose();
	math::matrix(K, 12, 12).span(0, 3) = -A1 * (Ku * Sr - Fr) * A1.transpose();
	math::matrix(K, 12, 12).span(6, 3) = +A1 * (Ku * Sr - Fr) * A1.transpose();
	math::matrix(K, 12, 12).span(3, 3) = +A1 * (Kt - Sr * Ku * Sr + Mr + Sr * Fr) * A1.transpose();
	math::matrix(K, 12, 12).span(0, 3) = math::matrix(K, 12, 12).span(0, 3) * t1.rotation_gradient();
	math::matrix(K, 12, 12).span(3, 3) = math::matrix(K, 12, 12).span(3, 3) * t1.rotation_gradient();
	math::matrix(K, 12, 12).span(6, 3) = math::matrix(K, 12, 12).span(6, 3) * t1.rotation_gradient();
	math::matrix(K, 12, 12).span(9, 3) = math::matrix(K, 12, 12).span(9, 3) * t1.rotation_gradient();
	math::matrix(K, 12, 12).span(3, 9) = math::matrix(K, 12, 12).span(3, 9) * t2.rotation_gradient();
	math::matrix(K, 12, 12).span(9, 9) = math::matrix(K, 12, 12).span(9, 9) * t2.rotation_gradient();
}
static void test_force(void)
{
	double d[12];
	math::vector fia(12);
	math::vector fin(12);
	math::vector fie(12);
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		math::vector(d, 12).randu();
		internal_force_2(fia.data(), d, nullptr);
		math::ndiff(energy, fin.data(), d, nullptr, 1, 12, 1e-5);
		fie = fia - fin;
		if(fie.norm() < 1e-5)
		{
			printf("%04d: ok!\n", i);
		}
		else
		{
			fie.print("fie", 1e-5);
			break;
		}
	}
}
static void test_stiffness(void)
{
	double d[12];
	math::matrix Ka(12, 12);
	math::matrix Kn(12, 12);
	math::matrix Ke(12, 12);
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		math::vector(d, 12).randu();
		stiffness(Ka.data(), d, nullptr);
		math::ndiff(internal_force_1, Kn.data(), d, nullptr, 12, 12, 1e-5);
		Ke = Ka - Kn;
		if(Ke.norm() < 1e-5)
		{
			printf("%04d: ok!\n", i);
		}
		else
		{
			Ke.print("Ke", 1e-5);
			break;
		}
	}
}

void tests::fem::revolute(void)
{
	true ? test_force() : test_stiffness();
}
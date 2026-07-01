//std
#include <ctime>
#include <cstdio>
#include <cstdlib>

//Math
#include "Math/inc/Miscellaneous/util.hpp"
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Mat3.hpp"
#include "Math/inc/Linear/Quat.hpp"

//Test
#include "Math/Test/inc/fem.hpp"

static double Lr;
static math::Vec3 z1, z2;
static math::Quat q1r, q2r;
static math::Vec3 u1, u2, u3;

static uint32_t nt = 10000;
static uint32_t what, order;
static math::Vector d(12), ga(12), gn(12), ge(12);
static math::Matrix ha(12, 12), hn(12, 12), he(12, 12);

/*
	(x) Rotation
	R1 = exp(t1) * R1n
	R2 = exp(t2) * R2n
	(x) Triads
	ak = R1 * uk
	bk = R2 * uk
	(x) Constraints
	f1 = a1 . b3 = 0
	f2 = a2 . b3 = 0
	d1 = a1 . (x2 - x1) = 0
	d2 = a2 . (x2 - x1) = 0
	d3 = a3 . (x2 - x1) - Lr = 0
*/

static void setup(void)
{
	z1.randu();
	z2.randu();
	q1r.randu();
	q2r.randu();
	u3 = (z2 - z1).unit();
	Lr = (z2 - z1).norm();
	u3.triad(u1, u2);
}

static void constraint_1(double* constraint, const double* d, void** args)
{
	//data
	const math::Quat q1 = math::Vec3(d + 3).quaternion() * q1r;
	const math::Quat q2 = math::Vec3(d + 9).quaternion() * q2r;
	//triads
	const math::Vec3 a1 = q1.rotate(u1);
	const math::Vec3 b3 = q2.rotate(u3);
	//constraint
	constraint[0] = a1.inner(b3);
}
static void constraint_2(double* constraint, const double* d, void** args)
{
	//data
	const math::Vec3 x1 = z1 + math::Vec3(d + 0);
	const math::Vec3 x2 = z1 + math::Vec3(d + 6);
	const math::Quat q1 = math::Vec3(d + 3).quaternion() * q1r;
	const math::Quat q2 = math::Vec3(d + 9).quaternion() * q2r;
	//triads
	const math::Vec3 a1 = q1.rotate(u1);
	const math::Vec3 a2 = q1.rotate(u2);
	const math::Vec3 a3 = q1.rotate(u3);
	const math::Vec3 b3 = q2.rotate(u3);
	//constraint
	constraint[0] = a2.inner(b3);
}
static void constraint_3(double* constraint, const double* d, void** args)
{
	//data
	const math::Vec3 x1 = z1 + math::Vec3(d + 0);
	const math::Vec3 x2 = z2 + math::Vec3(d + 6);
	const math::Quat q1 = math::Vec3(d + 3).quaternion() * q1r;
	//triads
	const math::Vec3 a1 = q1.rotate(u1);
	//constraint
	constraint[0] = a1.inner(x2 - x1);
}
static void constraint_4(double* constraint, const double* d, void** args)
{
	//data
	const math::Vec3 x1 = z1 + math::Vec3(d + 0);
	const math::Vec3 x2 = z2 + math::Vec3(d + 6);
	const math::Quat q1 = math::Vec3(d + 3).quaternion() * q1r;
	//triads
	const math::Vec3 a2 = q1.rotate(u2);
	//constraints
	constraint[0] = a2.inner(x2 - x1);
}
static void constraint_5(double* constraint, const double* d, void** args)
{
	//data
	const math::Vec3 x1 = z1 + math::Vec3(d + 0);
	const math::Vec3 x2 = z2 + math::Vec3(d + 6);
	const math::Quat q1 = math::Vec3(d + 3).quaternion() * q1r;
	//triads
	const math::Vec3 a3 = q1.rotate(u3);
	//constraint
	constraint[0] = a3.inner(x2 - x1) - Lr;
}
void constraint(double* constraint, const double* d, void** args)
{
	//data
	void(*constraints[])(double*, const double*, void**) = {
		constraint_1, constraint_2, constraint_3, constraint_4, constraint_5
	};
	//constraint
	constraints[what](constraint, d, args);
}

static void gradient_1(double* gradient, const double* d, void** args)
{
	//data
	const math::Vec3 t1(d + 3);
	const math::Vec3 t2(d + 9);
	const math::Quat q1 = t1.quaternion() * q1r;
	const math::Quat q2 = t2.quaternion() * q2r;
	//triads
	const math::Vec3 a1 = q1.rotate(u1);
	const math::Vec3 b3 = q2.rotate(u3);
	//gradient
	math::Vector(gradient, 12).zeros();
	math::Vector(gradient, 12).span(3, 0, 3, 1) = t1.rotation_gradient(a1.cross(b3), true);
	math::Vector(gradient, 12).span(9, 0, 3, 1) = t2.rotation_gradient(b3.cross(a1), true);
}
static void gradient_2(double* gradient, const double* d, void** args)
{
	//data
	const math::Vec3 t1(d + 3);
	const math::Vec3 t2(d + 9);
	const math::Quat q1 = t1.quaternion() * q1r;
	const math::Quat q2 = t2.quaternion() * q2r;
	//triads
	const math::Vec3 a2 = q1.rotate(u2);
	const math::Vec3 b3 = q2.rotate(u3);
	//gradient
	math::Vector(gradient, 12).zeros();
	math::Vector(gradient, 12).span(3, 0, 3, 1) = t1.rotation_gradient(a2.cross(b3), true);
	math::Vector(gradient, 12).span(9, 0, 3, 1) = t2.rotation_gradient(b3.cross(a2), true);
}
static void gradient_3(double* gradient, const double* d, void** args)
{
	//data
	const math::Vec3 t1(d + 3);
	const math::Quat q1 = t1.quaternion() * q1r;
	const math::Vec3 x1 = z1 + math::Vec3(d + 0);
	const math::Vec3 x2 = z2 + math::Vec3(d + 6);
	//triads
	const math::Vec3 a1 = q1.rotate(u1);
	//gradient
	math::Vector(gradient, 12).zeros();
	math::Vector(gradient, 12).span(0, 0, 3, 1) = -a1;
	math::Vector(gradient, 12).span(6, 0, 3, 1) = +a1;
	math::Vector(gradient, 12).span(3, 0, 3, 1) = t1.rotation_gradient(a1.cross(x2 - x1), true);
}
static void gradient_4(double* gradient, const double* d, void** args)
{
	//data
	const math::Vec3 t1(d + 3);
	const math::Quat q1 = t1.quaternion() * q1r;
	const math::Vec3 x1 = z1 + math::Vec3(d + 0);
	const math::Vec3 x2 = z2 + math::Vec3(d + 6);
	//triads
	const math::Vec3 a2 = q1.rotate(u2);
	//gradient
	math::Vector(gradient, 12).zeros();
	math::Vector(gradient, 12).span(0, 0, 3, 1) = -a2;
	math::Vector(gradient, 12).span(6, 0, 3, 1) = +a2;
	math::Vector(gradient, 12).span(3, 0, 3, 1) = t1.rotation_gradient(a2.cross(x2 - x1), true);
}
static void gradient_5(double* gradient, const double* d, void** args)
{
	//data
	const math::Vec3 t1(d + 3);
	const math::Quat q1 = t1.quaternion() * q1r;
	const math::Vec3 x1 = z1 + math::Vec3(d + 0);
	const math::Vec3 x2 = z2 + math::Vec3(d + 6);
	//triads
	const math::Vec3 a3 = q1.rotate(u3);
	//gradient
	math::Vector(gradient, 12).zeros();
	math::Vector(gradient, 12).span(0, 0, 3, 1) = -a3;
	math::Vector(gradient, 12).span(6, 0, 3, 1) = +a3;
	math::Vector(gradient, 12).span(3, 0, 3, 1) = t1.rotation_gradient(a3.cross(x2 - x1), true);
}
void gradient(double* gradient, const double* d, void** args)
{
	//data
	void(*gradients[])(double*, const double*, void**) = {
		gradient_1, gradient_2, gradient_3, gradient_4, gradient_5
	};
	//gradient
	gradients[what](gradient, d, args);
}

static void hessian_1(double* hessian, const double* d, void** args)
{
	//data
	const math::Vec3 t1(d + 3);
	const math::Vec3 t2(d + 9);
	const math::Quat q1 = t1.quaternion() * q1r;
	const math::Quat q2 = t2.quaternion() * q2r;
	const math::Mat3 T1 = t1.rotation_gradient();
	const math::Mat3 T2 = t2.rotation_gradient();
	//triads
	const math::Vec3 a1 = q1.rotate(u1);
	const math::Vec3 b3 = q2.rotate(u3);
	const math::Mat3 A11 = t1.rotation_hessian(a1.cross(b3), true);
	const math::Mat3 A21 = t2.rotation_hessian(b3.cross(a1), true);
	//hessian
	math::Matrix(hessian, 12, 12).zeros();
	math::Matrix(hessian, 12, 12).span(3, 9) = -T1.transpose() * a1.spin() * b3.spin() * T2;
	math::Matrix(hessian, 12, 12).span(9, 3) = -T2.transpose() * b3.spin() * a1.spin() * T1;
	math::Matrix(hessian, 12, 12).span(3, 3) = +T1.transpose() * b3.spin() * a1.spin() * T1 + A11;
	math::Matrix(hessian, 12, 12).span(9, 9) = +T2.transpose() * a1.spin() * b3.spin() * T2 + A21;
}
static void hessian_2(double* hessian, const double* d, void** args)
{
	//data
	const math::Vec3 t1(d + 3);
	const math::Vec3 t2(d + 9);
	const math::Quat q1 = t1.quaternion() * q1r;
	const math::Quat q2 = t2.quaternion() * q2r;
	const math::Mat3 T1 = t1.rotation_gradient();
	const math::Mat3 T2 = t2.rotation_gradient();
	//triads
	const math::Vec3 a2 = q1.rotate(u2);
	const math::Vec3 b3 = q2.rotate(u3);
	const math::Mat3 A11 = t1.rotation_hessian(a2.cross(b3), true);
	const math::Mat3 A21 = t2.rotation_hessian(b3.cross(a2), true);
	//hessian
	math::Matrix(hessian, 12, 12).zeros();
	math::Matrix(hessian, 12, 12).span(3, 9) = -T1.transpose() * a2.spin() * b3.spin() * T2;
	math::Matrix(hessian, 12, 12).span(9, 3) = -T2.transpose() * b3.spin() * a2.spin() * T1;
	math::Matrix(hessian, 12, 12).span(3, 3) = +T1.transpose() * b3.spin() * a2.spin() * T1 + A11;
	math::Matrix(hessian, 12, 12).span(9, 9) = +T2.transpose() * a2.spin() * b3.spin() * T2 + A21;
}
static void hessian_3(double* hessian, const double* d, void** args)
{
	//data
	const math::Vec3 t1(d + 3);
	const math::Quat q1 = t1.quaternion() * q1r;
	const math::Mat3 T1 = t1.rotation_gradient();
	const math::Vec3 x1 = z1 + math::Vec3(d + 0);
	const math::Vec3 x2 = z2 + math::Vec3(d + 6);
	//triads
	const math::Vec3 a1 = q1.rotate(u1);
	const math::Mat3 A13 = t1.rotation_hessian(a1.cross(x2 - x1), true);
	//hessian
	math::Matrix(hessian, 12, 12).zeros();
	math::Matrix(hessian, 12, 12).span(0, 3) = +a1.spin() * T1;
	math::Matrix(hessian, 12, 12).span(6, 3) = -a1.spin() * T1;
	math::Matrix(hessian, 12, 12).span(3, 0) = -T1.transpose() * a1.spin();
	math::Matrix(hessian, 12, 12).span(3, 6) = +T1.transpose() * a1.spin();
	math::Matrix(hessian, 12, 12).span(3, 3) = +T1.transpose() * (x2 - x1).spin() * a1.spin() * T1 + A13;
}
static void hessian_4(double* hessian, const double* d, void** args)
{
	//data
	const math::Vec3 t1(d + 3);
	const math::Quat q1 = t1.quaternion() * q1r;
	const math::Mat3 T1 = t1.rotation_gradient();
	const math::Vec3 x1 = z1 + math::Vec3(d + 0);
	const math::Vec3 x2 = z2 + math::Vec3(d + 6);
	//triads
	const math::Vec3 a2 = q1.rotate(u2);
	const math::Mat3 A14 = t1.rotation_hessian(a2.cross(x2 - x1), true);
	//hessian
	math::Matrix(hessian, 12, 12).zeros();
	math::Matrix(hessian, 12, 12).span(0, 3) = +a2.spin() * T1;
	math::Matrix(hessian, 12, 12).span(6, 3) = -a2.spin() * T1;
	math::Matrix(hessian, 12, 12).span(3, 0) = -T1.transpose() * a2.spin();
	math::Matrix(hessian, 12, 12).span(3, 6) = +T1.transpose() * a2.spin();
	math::Matrix(hessian, 12, 12).span(3, 3) = +T1.transpose() * (x2 - x1).spin() * a2.spin() * T1 + A14;
}
static void hessian_5(double* hessian, const double* d, void** args)
{
	//data
	const math::Vec3 t1(d + 3);
	const math::Quat q1 = t1.quaternion() * q1r;
	const math::Mat3 T1 = t1.rotation_gradient();
	const math::Vec3 x1 = z1 + math::Vec3(d + 0);
	const math::Vec3 x2 = z2 + math::Vec3(d + 6);
	//triads
	const math::Vec3 a3 = q1.rotate(u3);
	const math::Mat3 A15 = t1.rotation_hessian(a3.cross(x2 - x1), true);
	//hessian
	math::Matrix(hessian, 12, 12).zeros();
	math::Matrix(hessian, 12, 12).span(0, 3) = +a3.spin() * T1;
	math::Matrix(hessian, 12, 12).span(6, 3) = -a3.spin() * T1;
	math::Matrix(hessian, 12, 12).span(3, 0) = -T1.transpose() * a3.spin();
	math::Matrix(hessian, 12, 12).span(3, 6) = +T1.transpose() * a3.spin();
	math::Matrix(hessian, 12, 12).span(3, 3) = +T1.transpose() * (x2 - x1).spin() * a3.spin() * T1 + A15;
}
static void hessian(double* hessian, const double* d, void** args)
{
	//data
	void(*hessians[])(double*, const double*, void**) = {
		hessian_1, hessian_2, hessian_3, hessian_4, hessian_5
	};
	//hessian
	hessians[what](hessian, d, args);
}

static void menu_what(void)
{
	while(true)
	{
		printf("What:\n");
		printf("(1) rotation 1\n");
		printf("(2) rotation 2\n");
		printf("(3) position 1\n");
		printf("(4) position 2\n");
		printf("(5) position 3\n");
		const int args = scanf("%d", &what);
		if(args == 1 && (what >= 1 || what <= 5)) break;
		printf("Invalid option!\n");
	}
	what -= 1;
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

static void test_gradient(void)
{
	srand((uint32_t) time(nullptr));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		gradient(ga.data(), d.data(), nullptr);
		math::ndiff(constraint, gn.data(), d.data(), nullptr, 1, 12, 1.00e-5);
		ge = ga - gn;
		const bool test = ge.norm() < 1e-5;
		printf("Test %04d: %s\n", i, test ? "ok" : "not ok");
		if(!test)
		{
			ga.print("ga");
			gn.print("gn");
			ge.print("ge", 1e-5);
			break;
		}
	}
}
static void test_hessian(void)
{
	srand((uint32_t) time(nullptr));
	for(uint32_t i = 0; i < nt; i++)
	{
		setup();
		hessian(ha.data(), d.data(), nullptr);
		math::ndiff(gradient, hn.data(), d.data(), nullptr, 12, 12, 1.00e-5);
		he = ha - hn;
		const bool test = he.norm() < 1e-5;
		printf("Test %04d: %s\n", i, test ? "ok" : "not ok");
		if(!test)
		{
			ha.print("ha");
			hn.print("hn");
			he.print("he", 1e-5);
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
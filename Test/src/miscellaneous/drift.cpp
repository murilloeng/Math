//std
#include <ctime>
#include <cstdio>
#include <cstdlib>

//Math
#include "Math/inc/Miscellaneous/util.hpp"
#include "Math/inc/Linear/Mat3.hpp"
#include "Math/inc/Linear/Vec3.hpp"

//Test
#include "Math/Test/inc/miscellaneous.hpp"

static double h = 0.01;
static const double g = 0.50;
static const double b = 0.25;

static math::Mat3 J;
static math::Vec3 qn, wn, an, t, dt, ddt;

static void fun(double* r, const double* x, void** args)
{
	math::Vec3 rm(r);
	const math::Vec3 t(x);
	const math::Vec3 ddt = an + (t - h * wn - h * h / 2 * an) / (h * h * b);
	const math::Vec3 dt = wn + h * an + g / (h * b) * (t - h * wn + h * h / 2 * an);
	const math::Vec3 w = t.rotation_gradient(dt, true);
	const math::Vec3 a = t.rotation_gradient(ddt, true) + t.rotation_hessian(dt, dt, true);
	rm = J * a + w.cross(J * w);
}
static void dfun(math::Mat3& Ka, const double* x, void** args)
{
	//data
	const math::Vec3 t(x);
	const math::Vec3 ddt = an + (t - h * wn - h * h / 2 * an) / (h * h * b);
	const math::Vec3 dt = wn + h * an + g / (h * b) * (t - h * wn + h * h / 2 * an);
	//velocity
	const math::Vec3 w = t.rotation_gradient(dt, true);
	const math::Mat3 Dw_Ddt = t.rotation_gradient(true);
	const math::Mat3 Dw_Dt = t.rotation_hessian(dt, true);
	const math::Mat3 Dw = Dw_Dt + g / (h * b) * Dw_Ddt;
	//acceleration
	const math::Mat3 Da_Dddt = t.rotation_gradient(true);
	const math::Mat3 Da_Dt = t.rotation_hessian(ddt, true) + t.rotation_third(dt, dt, true, true);
	const math::Mat3 Da_Ddt = t.rotation_hessian(dt, true) + t.rotation_third(dt, dt, true, false);
	const math::Mat3 Da = Da_Dt + g / (h * b) * Da_Ddt + 1 / (h * h * b) * Da_Dddt;
	//system
	Ka = J * Da + (w.spin() * J - (J * w).spin()) * Dw;
}

void tests::miscellaneous::drift(void)
{
	//test
	J[0] = 1;
	J[1] = 2;
	J[2] = 3;
	math::Mat3 Ka, Kn;
	const uint32_t n = 10000;
	srand((uint32_t) time(nullptr));
	for(uint32_t i = 0; i < n; i++)
	{
		t.randu();
		wn.randu();
		an.randu();
		dt.randu();
		ddt.randu();
		dfun(Ka, t.data(), nullptr);
		math::ndiff(fun, Kn.data(), t.data(), nullptr, 3, 3, 1e-8);
		const double r = (Ka - Kn).norm();
		if (r < 1e-2)
		{
			printf("%04d %+.2e ok\n", i, r);
		}
		else
		{
			printf("%04d %+.2e not ok\n", i, r);
			break;
		}
	}
}
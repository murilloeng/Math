//std
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>

//math
#include "Math/inc/misc/misc.hpp"
#include "Math/inc/linear/vec3.hpp"
#include "Math/inc/linear/mat3.hpp"

static double h = 0.01;
static const double g = 0.50;
static const double b = 0.25;

static math::mat3 J;
static math::vec3 qn, wn, an, t, dt, ddt;

void fun(double* r, const double* x)
{
	math::vec3 rm(r);
	const math::vec3 t(x);
	const math::vec3 ddt = an + (t - h * wn - h * h / 2 * an) / (h * h * b);
	const math::vec3 dt = wn + h * an + g / (h * b) * (t - h * wn + h * h / 2 * an);
	const math::vec3 w = t.rotation_gradient(dt, true);
	const math::vec3 a = t.rotation_gradient(ddt, true) + t.rotation_hessian(dt, dt, true);
	rm = J * a + w.cross(J * w);
}
void dfun(math::mat3& Ka, const double* x)
{
	//data
	const math::vec3 t(x);
	const math::vec3 ddt = an + (t - h * wn - h * h / 2 * an) / (h * h * b);
	const math::vec3 dt = wn + h * an + g / (h * b) * (t - h * wn + h * h / 2 * an);
	//velocity
	const math::vec3 w = t.rotation_gradient(dt, true);
	const math::mat3 Dw_Ddt = t.rotation_gradient(true);
	const math::mat3 Dw_Dt = t.rotation_hessian(dt, true);
	const math::mat3 Dw = Dw_Dt + g / (h * b) * Dw_Ddt;
	//acceleration
	const math::mat3 Da_Dddt = t.rotation_gradient(true);
	const math::mat3 Da_Dt = t.rotation_hessian(ddt, true) + t.rotation_higher(dt, dt, true, true);
	const math::mat3 Da_Ddt = t.rotation_hessian(dt, true) + t.rotation_higher(dt, dt, true, false);
	const math::mat3 Da = Da_Dt + g / (h * b) * Da_Ddt + 1 / (h * h * b) * Da_Dddt;
	//system
	Ka = J * Da + (w.spin() * J - (J * w).spin()) * Dw;
}

int main(void)
{
	//test
	J[0] = 1;
	J[1] = 2;
	J[2] = 3;
	math::mat3 Ka, Kn;
	srand(time(nullptr));
	const unsigned n = 10000;
	for(unsigned i = 0; i < n; i++)
	{
		t.randu();
		wn.randu();
		an.randu();
		dt.randu();
		ddt.randu();
		dfun(Ka, t.data());
		math::ndiff(fun, Kn.data(), t.data(), 3, 3, 1e-8);
		const double r = (Ka - Kn).norm();
		if(r < 1e-2)
		{
			printf("%04d %+.2e ok\n", i, r);
		}
		else
		{
			printf("%04d %+.2e not ok\n", i, r);
			break;
		}
	}
	//return
	return EXIT_SUCCESS;
}
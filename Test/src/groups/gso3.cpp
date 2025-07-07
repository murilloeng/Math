//std
#include <cmath>
#include <ctime>
#include <cstdio>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/groups/ASO3.hpp"
#include "Math/Math/inc/groups/GSO3.hpp"

//tests
#include "Math/Test/inc/tests.hpp"

//static
static void test_exponential(double* r, const double* v, const void** args)
{
	const math::vec3 a = (const double*) args[0];
	math::vec3(r + 0) = math::groups::ASO3(v).exponential() * a;
}
static void test_tangent(double* r, const double* v, const void** args)
{
	const math::vec3 a = (const double*) args[0];
	math::vec3(r + 0) = math::groups::ASO3(v).tangent(a);
}

//tests
void tests::groups::gso3_log(void)
{
	math::vec3 v, r;
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		v.randu();
		r = v - math::groups::ASO3(v).exponential().logarithm().vector();
		if(r.norm() > 1e-5)
		{
			v.print("v");
			r.print("r");
			break;
		}
		else
		{
			printf("Test %5d: ok!\n", i);
		}
	}
}
void tests::groups::gso3_inverse(void)
{
	math::quat q, r;
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		q.randu();
		r = q * math::groups::GSO3(q).inverse().quaternion();
		r[0] -= 1;
		if(r.norm() > 1e-5)
		{
			q.print("v");
			r.print("r");
			break;
		}
		else
		{
			printf("Test %5d: ok!\n", i);
		}
	}
}
void tests::groups::aso3_tangent(void)
{
	math::vec3 a, v, r;
	const uint32_t nt = 100000;
	math::mat3 Ka, Kn, Kr, R, T;
	srand(uint32_t(time(nullptr)));
	const void* args[] = {a.data()};
	for(uint32_t i = 0; i < nt; i++)
	{
		a.randu();
		v.randu();
		T = math::groups::ASO3(v).tangent();
		R = math::groups::ASO3(v).exponential().matrix();
		Ka = -R * a.spin() * T;
		math::ndiff(test_exponential, Kn.data(), v.data(), args, 3, 3, 1e-5);
		Kr = Ka - Kn;
		if(Kr.norm() > 1e-5 * Ka.norm())
		{
			Ka.print("Ka");
			Kn.print("Kn");
			Kr.print("Kr");
			break;
		}
		else
		{
			printf("Test %5d: ok!\n", i);
		}
	}
}
void tests::groups::aso3_tangent_inverse(void)
{
	math::vec3 v;
	math::mat3 R, T, A;
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		v.randu();
		T = math::groups::ASO3(v).tangent();
		A = math::groups::ASO3(v).tangent_inverse();
		R = T * A - math::mat3::eye();
		if(R.norm() > 1e-5)
		{
			T.print("T");
			A.print("A");
			R.print("R");
			break;
		}
		else
		{
			printf("Test %5d: ok!\n", i);
		}
	}
}
void tests::groups::aso3_tangent_increment(void)
{
	math::vec3 a, v, r;
	math::mat3 Ka, Kn, Kr;
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	const void* args[] = {a.data()};
	for(uint32_t i = 0; i < nt; i++)
	{
		a.randu();
		v.randu();
		Ka = math::groups::ASO3(v).tangent_increment(a);
		math::ndiff(test_tangent, Kn.data(), v.data(), args, 3, 3, 1e-5);
		Kr = Ka - Kn;
		if(Kr.norm() > 1e-5 * Ka.norm())
		{
			Ka.print("Ka");
			Kn.print("Kn");
			Kr.print("Kr");
			break;
		}
		else
		{
			printf("Test %5d: ok!\n", i);
		}
	}
}
void tests::groups::aso3_tangent_inverse_increment(void)
{

}
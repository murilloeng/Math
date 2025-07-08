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
static void test_tangent_inverse(double* r, const double* v, const void** args)
{
	const math::vec3 a = (const double*) args[0];
	math::vec3(r + 0) = math::groups::ASO3(v).tangent_inverse(a);
}

//tests
void tests::groups::GSO3::log(void)
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
			printf("Test gso3 logarithm %5d: ok!\n", i);
		}
	}
}
void tests::groups::GSO3::inverse(void)
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
			printf("Test gso3 inverse %5d: ok!\n", i);
		}
	}
}
void tests::groups::GSO3::tangent(void)
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
		R = math::groups::ASO3(v).exponential();
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
			printf("Test GSO3 tangent %5d: ok!\n", i);
		}
	}
}
void tests::groups::GSO3::tangent_inverse(void)
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
			printf("Test GSO3 tangent inverse %5d: ok!\n", i);
		}
	}
}
void tests::groups::GSO3::tangent_indentity(void)
{
	math::mat3 Kr;
	math::groups::ASO3 object;
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		object.vector().randu();
		Kr = math::mat3(object.exponential()) * object.tangent() - object.tangent().transpose();
		if(Kr.norm() > 1e-5)
		{
			object.vector().print("v");
			Kr.print("Kr");
			break;
		}
		else
		{
			printf("Test GSO3 tangent identity %5d: ok!\n", i);
		}
	}
}
void tests::groups::GSO3::tangent_increment(void)
{
	math::vec3 a, v;
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
			printf("Test GSO3 tangent increment %5d: ok!\n", i);
		}
	}
}
void tests::groups::GSO3::tangent_inverse_increment(void)
{
	math::vec3 a, v;
	math::mat3 Ka, Kn, Kr;
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	const void* args[] = {a.data()};
	for(uint32_t i = 0; i < nt; i++)
	{
		a.randu();
		v.randu();
		Ka = math::groups::ASO3(v).tangent_inverse_increment(a);
		math::ndiff(test_tangent_inverse, Kn.data(), v.data(), args, 3, 3, 1e-5);
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
			printf("Test GSO3 tangent inverse increment %5d: ok!\n", i);
		}
	}
}
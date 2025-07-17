//std
#include <cmath>
#include <ctime>
#include <cstdio>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat3.hpp"
#include "Math/Math/inc/groups/ASE3.hpp"
#include "Math/Math/inc/groups/GSE3.hpp"

//tests
#include "Math/Test/inc/tests.hpp"

//static
static void test_exponential(double* r, const double* v, const void** args)
{
	r[3] = 1;
	const math::vec3 u(v + 0), w(v + 3);
	const math::vec3 a = (const double*) args[0];
	math::vec3(r + 0) = math::groups::ASE3(u, w).exponential() * a;
}

//tests
void tests::groups::GSE3::log(void)
{
	math::vec3 u, w, ru, rw;
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		u.randu();
		w.randu();
		math::groups::ASE3 object(u, w);
		ru = u - object.exponential().logarithm().vector_u();
		rw = w - object.exponential().logarithm().vector_w();
		if(ru.norm() > 1e-5 || rw.norm() > 1e-5)
		{
			u.print("u");
			w.print("w");
			ru.print("ru");
			rw.print("rw");
			break;
		}
		else
		{
			printf("Test GSE3 logarithm %5d: ok!\n", i);
		}
	}
}
void tests::groups::GSE3::inverse(void)
{
	math::groups::GSE3 g, r;
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		g.vector().randu();
		g.quaternion().randu();
		r = g * g.inverse();
		r.quaternion()[0] -= 1;
		if(r.vector().norm() > 1e-5 || r.quaternion().norm() > 1e-5)
		{
			r.vector().print("rv");
			r.quaternion().print("rq");
			break;
		}
		else
		{
			printf("Test GSE3 inverse %5d: ok!\n", i);
		}
	}
}
void tests::groups::GSE3::tangent(void)
{
	const uint32_t nt = 100000;
	math::vector a(3), v(6), r(4);
	srand(uint32_t(time(nullptr)));
	const void* args[] = {a.data()};
	math::matrix Ka(4, 6), Kn(4, 6), Kr(4, 6), A(4, 6, math::mode::zeros);
	for(uint32_t i = 0; i < nt; i++)
	{
		a.randu();
		v.randu();
		const math::vec3 u = v.data() + 0;
		const math::vec3 w = v.data() + 3;
		A.span(0, 0, 3, 3) = math::matrix::eye(3, 3);
		// A.span(0, 3, 3, 3) = -((math::matrix&) math::vec3(a.data()).spin());
		const math::matrix T = math::groups::ASE3(u, w).tangent();
		const math::mat4 H = math::groups::ASE3(u, w).exponential();
		Ka = ((math::matrix&) H) * A * T;
		math::ndiff(test_exponential, Kn.data(), v.data(), args, 4, 6, 1e-5);
		Kr = Ka - Kn;
		if(Kr.norm() > 1e-5 * Ka.norm())
		{
			Ka.print("Ka");
			Kn.print("Kn");
			Kr.print("Kr", 1e-5);
			break;
		}
		else
		{
			printf("Test GSE3 tangent %5d: ok!\n", i);
		}
	}
}
void tests::groups::GSE3::tangent_inverse(void)
{
	math::groups::ASE3 object;
	const uint32_t nt = 100000;
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		object.vector_u().randu();
		object.vector_w().randu();
		const math::matrix T = object.tangent();
		const math::matrix Ti = object.tangent_inverse();
		const math::matrix Tr = Ti * T - math::matrix::eye(6, 6);
		if(Tr.norm() > 1e-5)
		{
			T.print("T");
			Ti.print("Ti");
			Tr.print("Tr");
			break;
		}
		else
		{
			printf("Test GSE3 tangent inverse %5d: ok!\n", i);
		}
	}
}
void tests::groups::GSE3::tangent_increment(void)
{
	// math::vec3 a, v;
	// math::mat3 Ka, Kn, Kr;
	// const uint32_t nt = 100000;
	// srand(uint32_t(time(nullptr)));
	// const void* args[] = {a.data()};
	// for(uint32_t i = 0; i < nt; i++)
	// {
	// 	a.randu();
	// 	v.randu();
	// 	Ka = math::groups::ASO3(v).tangent_increment(a);
	// 	math::ndiff(test_tangent, Kn.data(), v.data(), args, 3, 3, 1e-5);
	// 	Kr = Ka - Kn;
	// 	if(Kr.norm() > 1e-5 * Ka.norm())
	// 	{
	// 		Ka.print("Ka");
	// 		Kn.print("Kn");
	// 		Kr.print("Kr");
	// 		break;
	// 	}
	// 	else
	// 	{
	// 		printf("Test GSE3 tangent increment %5d: ok!\n", i);
	// 	}
	// }
}
void tests::groups::GSE3::tangent_inverse_increment(void)
{
	// math::vec3 a, v;
	// math::mat3 Ka, Kn, Kr;
	// const uint32_t nt = 100000;
	// srand(uint32_t(time(nullptr)));
	// const void* args[] = {a.data()};
	// for(uint32_t i = 0; i < nt; i++)
	// {
	// 	a.randu();
	// 	v.randu();
	// 	Ka = math::groups::ASO3(v).tangent_inverse_increment(a);
	// 	math::ndiff(test_tangent_inverse, Kn.data(), v.data(), args, 3, 3, 1e-5);
	// 	Kr = Ka - Kn;
	// 	if(Kr.norm() > 1e-5 * Ka.norm())
	// 	{
	// 		Ka.print("Ka");
	// 		Kn.print("Kn");
	// 		Kr.print("Kr");
	// 		break;
	// 	}
	// 	else
	// 	{
	// 		printf("Test GSE3 tangent inverse increment %5d: ok!\n", i);
	// 	}
	// }
}
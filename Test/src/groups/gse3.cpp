//std
#include <cmath>
#include <ctime>
#include <cstdio>

//Math
#include "Math/inc/Miscellaneous/util.hpp"
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Mat3.hpp"
#include "Math/inc/Groups/ASE3.hpp"
#include "Math/inc/Groups/GSE3.hpp"

//tests
#include "Math/Test/inc/groups.hpp"

//static
static void test_exponential(double* r, const double* v, const void** args)
{
	r[3] = 1;
	const math::Vec3 u(v + 0), w(v + 3);
	const math::Vec3 a = (const double*) args[0];
	math::Vec3(r + 0) = math::groups::ASE3(u, w).exponential() * a;
}

//tests
void tests::groups::GSE3::log(void)
{
	math::Vec3 u, w, ru, rw;
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
		g.Vector().randu();
		g.quaternion().randu();
		r = g * g.inverse();
		r.quaternion()[0] -= 1;
		if(r.Vector().norm() > 1e-5 || r.quaternion().norm() > 1e-5)
		{
			r.Vector().print("rv");
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
	math::Vector a(3), v(6), r(4);
	srand(uint32_t(time(nullptr)));
	const void* args[] = {a.data()};
	math::Matrix Ka(4, 6), Kn(4, 6), Kr(4, 6), A(4, 6, math::mode::zeros);
	for(uint32_t i = 0; i < nt; i++)
	{
		a.randu();
		v.randu();
		const math::Vec3 u = v.data() + 0;
		const math::Vec3 w = v.data() + 3;
		A.Span(0, 0, 3, 3) = math::Matrix::eye(3, 3);
		// A.Span(0, 3, 3, 3) = -((math::Matrix&) math::Vec3(a.data()).spin());
		const math::Matrix T = math::groups::ASE3(u, w).tangent();
		const math::Mat4 H = math::groups::ASE3(u, w).exponential();
		Ka = ((math::Matrix&) H) * A * T;
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
		const math::Matrix T = object.tangent();
		const math::Matrix Ti = object.tangent_inverse();
		const math::Matrix Tr = Ti * T - math::Matrix::eye(6, 6);
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
	// math::Vec3 a, v;
	// math::Mat3 Ka, Kn, Kr;
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
	// math::Vec3 a, v;
	// math::Mat3 Ka, Kn, Kr;
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
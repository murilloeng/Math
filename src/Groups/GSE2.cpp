//std
#include <cmath>
#include <cstdio>

//Math
#include "Math/inc/Groups/ASE2.hpp"
#include "Math/inc/Groups/GSE2.hpp"
#include "Math/inc/Miscellaneous/rotations.hpp"

namespace math
{
	//constructor
	GSE2::GSE2(const double* data) : m_data{data[0], data[1], data[2]}
	{
		return;
	}
	GSE2::GSE2(double x1, double x2, double t) : m_data{x1, x2, t}
	{
		return;
	}

	//random
	GSE2& GSE2::randu(double a, double b)
	{
		return *this = ASE2().randu(a, b).exp();
	}

	//actions
	Matrix GSE2::matrix(void) const
	{
		const double x1 = m_data[0];
		const double x2 = m_data[1];
		const double ct = cos(m_data[2]);
		const double st = sin(m_data[2]);
		return Matrix({{ct, -st, x1}, {st, ct, x2}, {0, 0, 1}});
	}
	Vector GSE2::action(const Vector& v) const
	{
		const double x1 = m_data[0];
		const double x2 = m_data[1];
		const double ct = cos(m_data[2]);
		const double st = sin(m_data[2]);
		return Vector({x1 + ct * v[0] - st * v[1], x2 + st * v[0] + ct * v[1]});
	}

	//operations
	ASE2 GSE2::log(void) const
	{
		const double x1 = m_data[0];
		const double x2 = m_data[1];
		const double tr = std::remainder(m_data[2], 2 * M_PI);
		const double fr = fn(tr, 1) / fn(tr, 2) / 2;
		return ASE2(fr * x1 + tr * x2 / 2, fr * x2 - tr * x1 / 2, tr);
	}
	GSE2 GSE2::inverse(void) const
	{
		const double x1 = m_data[0];
		const double x2 = m_data[1];
		const double ct = cos(m_data[2]);
		const double st = sin(m_data[2]);
		return GSE2(-ct * x1 - st * x2, st * x1 - ct * x2, -m_data[2]);
	}
	GSE2 GSE2::operator*(const GSE2& element) const
	{
		const double x1 = m_data[0];
		const double x2 = m_data[1];
		const double ct = cos(m_data[2]);
		const double st = sin(m_data[2]);
		const double z1 = element.m_data[0];
		const double z2 = element.m_data[1];
		return GSE2(x1 + ct * z1 - st * z2, x2 + st * z1 + ct * z2, m_data[2] + element.m_data[2]);
	}
	GSE2 GSE2::interpolate(const GSE2& element, double t) const
	{
		return *this * (t * (inverse() * element).log()).exp();
	}

	//tests
	bool GSE2::test_matrix(uint32_t nt, double tolerance)
	{
		//test
		GSE2 g1, g2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			g1.randu(-1, +1);
			g2.randu(-1, +1);
			test = test && ((g1 * g2).matrix() - g1.matrix() * g2.matrix()).norm() < tolerance;
		}
		//print
		printf("GSE2 matrix test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool GSE2::test_inverse(uint32_t nt, double tolerance)
	{
		//test
		GSE2 g;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			g.randu(-1, +1);
			test = test && (g * g.inverse()).log().vector().norm() < tolerance;
		}
		//print
		printf("GSE2 inverse test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool GSE2::test_log_exp(uint32_t nt, double tolerance)
	{
		//test
		GSE2 g;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			g.randu(-1, +1);
			test = test && (g.matrix() - g.log().exp().matrix()).norm() < tolerance;
		}
		//print
		printf("GSE2 log-exp test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
}
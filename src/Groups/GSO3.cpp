//std
#include <cmath>
#include <cstdio>

//Math
#include "Math/inc/Groups/ASO3.hpp"
#include "Math/inc/Groups/GSO3.hpp"
#include "Math/inc/Miscellaneous/rotations.hpp"

namespace math
{
	//constructor
	GSO3::GSO3(double t1, double t2, double t3)
	{
		//data
		const double tn = sqrt(t1 * t1 + t2 * t2 + t3 * t3);
		//quaternion
		m_data[0] = fn(tn / 2, 0);
		m_data[1] = fn(tn / 2, 1) * t1 / 2;
		m_data[2] = fn(tn / 2, 1) * t2 / 2;
		m_data[3] = fn(tn / 2, 1) * t3 / 2;
	}
	GSO3::GSO3(double q0, double q1, double q2, double q3) : m_data{q0, q1, q2, q3}
	{
		return;
	}

	//random
	GSO3& GSO3::randu(double a, double b)
	{
		return *this = ASO3().randu(a, b).exp();
	}

	//actions
	Matrix GSO3::matrix(void) const
	{
		//data
		Matrix R(3, 3);
		const double q0 = m_data[0];
		const double q1 = m_data[1];
		const double q2 = m_data[2];
		const double q3 = m_data[3];
		//matrix
		R[1 + 3 * 0] = 2 * (q1 * q2 + q0 * q3);
		R[2 + 3 * 0] = 2 * (q1 * q3 - q0 * q2);
		R[0 + 3 * 1] = 2 * (q1 * q2 - q0 * q3);
		R[2 + 3 * 1] = 2 * (q2 * q3 + q0 * q1);
		R[0 + 3 * 2] = 2 * (q1 * q3 + q0 * q2);
		R[1 + 3 * 2] = 2 * (q2 * q3 - q0 * q1);
		R[0 + 3 * 0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
		R[1 + 3 * 1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
		R[2 + 3 * 2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
		//return
		return R;
	}
	Vector GSO3::action(const Vector& v) const
	{
		//data
		Vector r(3);
		const double v1 = v[0];
		const double v2 = v[1];
		const double v3 = v[2];
		const double q0 = m_data[0];
		const double q1 = m_data[1];
		const double q2 = m_data[2];
		const double q3 = m_data[3];
		const double b = 2 * (q1 * v1 + q2 * v2 + q3 * v3);
		const double a = q0 * q0 - q1 * q1 - q2 * q2 - q3 * q3;
		//vector
		r[0] = a * v1 + b * q1 + 2 * q0 * (q2 * v3 - q3 * v2);
		r[1] = a * v2 + b * q2 + 2 * q0 * (q3 * v1 - q1 * v3);
		r[2] = a * v3 + b * q3 + 2 * q0 * (q1 * v2 - q2 * v1);
		//return
		return r;
	}

	//operations
	ASO3 GSO3::log(void) const
	{
		//data
		const double q0 = m_data[0];
		const double q1 = m_data[1];
		const double q2 = m_data[2];
		const double q3 = m_data[3];
		const double tn = atan2(sqrt(q1 * q1 + q2 * q2 + q3 * q3), q0);
		//return
		const double f1 = fn(tn, 1);
		return ASO3(2 * q1 / f1, 2 * q2 / f1, 2 * q3 / f1);
	}
	GSO3 GSO3::inverse(void) const
	{
		const double q0 = m_data[0];
		const double q1 = m_data[1];
		const double q2 = m_data[2];
		const double q3 = m_data[3];
		return GSO3(q0, -q1, -q2, -q3);
	}
	GSO3 GSO3::operator*(const GSO3& element) const
	{
		//data
		GSO3 r;
		r.m_data[0] = m_data[0] * element.m_data[0] - m_data[1] * element.m_data[1] - m_data[2] * element.m_data[2] - m_data[3] * element.m_data[3];
		r.m_data[1] = m_data[0] * element.m_data[1] + m_data[1] * element.m_data[0] + m_data[2] * element.m_data[3] - m_data[3] * element.m_data[2];
		r.m_data[2] = m_data[0] * element.m_data[2] + m_data[2] * element.m_data[0] + m_data[3] * element.m_data[1] - m_data[1] * element.m_data[3];
		r.m_data[3] = m_data[0] * element.m_data[3] + m_data[3] * element.m_data[0] + m_data[1] * element.m_data[2] - m_data[2] * element.m_data[1];
		//return
		return r;
	}
	GSO3 GSO3::interpolate(const GSO3& element, double t) const
	{
		return *this * (t * (inverse() * element).log()).exp();
	}

	//tests
	bool GSO3::test_matrix(uint32_t nt, double tolerance)
	{
		//test
		GSO3 g1, g2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			g1.randu(-1, +1);
			g2.randu(-1, +1);
			test = test && ((g1 * g2).matrix() - g1.matrix() * g2.matrix()).norm() < tolerance;
		}
		//print
		printf("GSO3 matrix test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool GSO3::test_inverse(uint32_t nt, double tolerance)
	{
		//test
		GSO3 g;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			g.randu(-1, +1);
			test = test && (g * g.inverse()).log().vector().norm() < tolerance;
		}
		//print
		printf("GSO3 inverse test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool GSO3::test_log_exp(uint32_t nt, double tolerance)
	{
		//test
		GSO3 g;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			g.randu(-1, +1);
			test = test && (g.matrix() - g.log().exp().matrix()).norm() < tolerance;
		}
		//print
		printf("GSO3 log-exp test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
}
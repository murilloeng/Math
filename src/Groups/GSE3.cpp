//std
#include <cstdio>

//Math
#include "Math/inc/Groups/ASO3.hpp"
#include "Math/inc/Groups/GSO3.hpp"
#include "Math/inc/Groups/ASE3.hpp"
#include "Math/inc/Groups/GSE3.hpp"
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Miscellaneous/rotations.hpp"

namespace math
{
	//constructor
	GSE3::GSE3(const double* data) : m_data{data[0], data[1], data[2], data[3], data[4], data[5], data[6]}
	{
		return;
	}
	GSE3::GSE3(double x1, double x2, double x3, double q0, double q1, double q2, double q3) : m_data{x1, x2, x3, q0, q1, q2, q3}
	{
		return;
	}

	//random
	GSE3& GSE3::randu(double a, double b)
	{
		return *this = ASE3().randu(a, b).exp();
	}

	//actions
	GSO3 GSE3::gso3(void) const
	{
		return GSO3(m_data[3], m_data[4], m_data[5], m_data[6]);
	}
	Matrix GSE3::matrix(void) const
	{
		//data
		Matrix H(4, 4);
		const double x1 = m_data[0];
		const double x2 = m_data[1];
		const double x3 = m_data[2];
		const double q0 = m_data[3];
		const double q1 = m_data[4];
		const double q2 = m_data[5];
		const double q3 = m_data[6];
		//matrix
		H(3, 3) = 1;
		H(0, 3) = x1;
		H(1, 3) = x2;
		H(2, 3) = x3;
		H(3, 0) = H(3, 1) = H(3, 2) = 0;
		H(1, 0) = 2 * (q1 * q2 + q0 * q3);
		H(2, 0) = 2 * (q1 * q3 - q0 * q2);
		H(0, 1) = 2 * (q1 * q2 - q0 * q3);
		H(2, 1) = 2 * (q2 * q3 + q0 * q1);
		H(0, 2) = 2 * (q1 * q3 + q0 * q2);
		H(1, 2) = 2 * (q2 * q3 - q0 * q1);
		H(0, 0) = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
		H(1, 1) = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
		H(2, 2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
		//return
		return H;
	}
	Vector GSE3::action(const Vector& v) const
	{
		//data
		Vector r(3);
		const double v1 = v[0];
		const double v2 = v[1];
		const double v3 = v[2];
		const double x1 = m_data[0];
		const double x2 = m_data[1];
		const double x3 = m_data[2];
		const double q0 = m_data[3];
		const double q1 = m_data[4];
		const double q2 = m_data[5];
		const double q3 = m_data[6];
		const double b = 2 * (q1 * v1 + q2 * v2 + q3 * v3);
		const double a = q0 * q0 - q1 * q1 - q2 * q2 - q3 * q3;
		//action
		r[0] = x1 + a * v1 + b * q1 + 2 * q0 * (q2 * v3 - q3 * v2);
		r[1] = x2 + a * v2 + b * q2 + 2 * q0 * (q3 * v1 - q1 * v3);
		r[2] = x3 + a * v3 + b * q3 + 2 * q0 * (q1 * v2 - q2 * v1);
		//return
		return r;
	}

	//operations
	double& GSE3::operator[](uint32_t index)
	{
		return m_data[index];
	}
	const double& GSE3::operator[](uint32_t index) const
	{
		return m_data[index];
	}

	//group
	ASE3 GSE3::log(void) const
	{
		//data
		const Vec3 x(m_data + 0);
		const ASO3 t = GSO3(m_data + 3).log();
		const Vector u = t.jacobian_inverse().transpose() * x;
		//return
		return ASE3(u[0], u[1], u[2], t[0], t[1], t[2]);
	}
	GSE3 GSE3::inverse(void) const
	{
		//data
		const double x1 = m_data[0];
		const double x2 = m_data[1];
		const double x3 = m_data[2];
		const double q0 = m_data[3];
		const double q1 = m_data[4];
		const double q2 = m_data[5];
		const double q3 = m_data[6];
		//position
		const double b = 2 * (q1 * x1 + q2 * x2 + q3 * x3);
		const double a = q0 * q0 - q1 * q1 - q2 * q2 - q3 * q3;
		const double z1 = a * x1 + b * q1 - 2 * q0 * (q2 * x3 - q3 * x2);
		const double z2 = a * x2 + b * q2 - 2 * q0 * (q3 * x1 - q1 * x3);
		const double z3 = a * x3 + b * q3 - 2 * q0 * (q1 * x2 - q2 * x1);
		//return
		return GSE3(-z1, -z2, -z3, q0, -q1, -q2, -q3);
	}
	GSE3 GSE3::operator*(const GSE3& g) const
	{
		//data
		GSE3 r;
		const double x1 = m_data[0];
		const double x2 = m_data[1];
		const double x3 = m_data[2];
		const double q0 = m_data[3];
		const double q1 = m_data[4];
		const double q2 = m_data[5];
		const double q3 = m_data[6];
		const double z1 = g.m_data[0];
		const double z2 = g.m_data[1];
		const double z3 = g.m_data[2];
		const double p0 = g.m_data[3];
		const double p1 = g.m_data[4];
		const double p2 = g.m_data[5];
		const double p3 = g.m_data[6];
		const double b = 2 * (q1 * z1 + q2 * z2 + q3 * z3);
		const double a = q0 * q0 - q1 * q1 - q2 * q2 - q3 * q3;
		//operation
		r.m_data[3] = q0 * p0 - q1 * p1 - q2 * p2 - q3 * p3;
		r.m_data[4] = q0 * p1 + q1 * p0 + q2 * p3 - q3 * p2;
		r.m_data[5] = q0 * p2 + q2 * p0 + q3 * p1 - q1 * p3;
		r.m_data[6] = q0 * p3 + q3 * p0 + q1 * p2 - q2 * p1;
		r.m_data[0] = x1 + a * z1 + b * q1 + 2 * q0 * (q2 * z3 - q3 * z2);
		r.m_data[1] = x2 + a * z2 + b * q2 + 2 * q0 * (q3 * z1 - q1 * z3);
		r.m_data[2] = x3 + a * z3 + b * q3 + 2 * q0 * (q1 * z2 - q2 * z1);
		//return
		return r;
	}
	GSE3 GSE3::interpolate(const GSE3& element, double t) const
	{
		return *this * (t * (inverse() * element).log()).exp();
	}

	//tests
	bool GSE3::test_matrix(uint32_t nt, double tolerance)
	{
		//test
		GSE3 g1, g2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			g1.randu(-1, +1);
			g2.randu(-1, +1);
			test = test && ((g1 * g2).matrix() - g1.matrix() * g2.matrix()).norm() < tolerance;
		}
		//print
		printf("GSE3 matrix test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool GSE3::test_inverse(uint32_t nt, double tolerance)
	{
		//test
		GSE3 g;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			g.randu(-1, +1);
			test = test && (g * g.inverse()).log().vector().norm() < tolerance;
		}
		//print
		printf("GSE3 inverse test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool GSE3::test_log_exp(uint32_t nt, double tolerance)
	{
		//test
		GSE3 g;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			g.randu(-1, +1);
			test = test && (g.matrix() - g.log().exp().matrix()).norm() < tolerance;
		}
		//print
		printf("GSE3 log-exp test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
}
//std
#include <cmath>
#include <cstdio>

//Math
#include "Math/inc/Groups/ASE3.hpp"
#include "Math/inc/Groups/GSE3.hpp"
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Mat3.hpp"
#include "Math/inc/Miscellaneous/rotations.hpp"

namespace math
{
	//constructor
	ASE3::ASE3(const double* data) : m_data{data[0], data[1], data[2], data[3], data[4], data[5]}
	{
		return;
	}
	ASE3::ASE3(double u1, double u2, double u3, double t1, double t2, double t3) : m_data{u1, u2, u3, t1, t2, t3}
	{
		return;
	}

	//group
	GSE3 ASE3::exp(void) const
	{
		//data
		const double u1 = m_data[0];
		const double u2 = m_data[1];
		const double u3 = m_data[2];
		const double t1 = m_data[3];
		const double t2 = m_data[4];
		const double t3 = m_data[5];
		const double tn = sqrt(t1 * t1 + t2 * t2 + t3 * t3);
		//functions
		const double f2 = fn(tn, 2);
		const double f3 = fn(tn, 3);
		const double f0 = fn(tn / 2, 0);
		const double f1 = fn(tn / 2, 1);
		//quaternion
		const double q0 = f0;
		const double q1 = f1 * t1 / 2;
		const double q2 = f1 * t2 / 2;
		const double q3 = f1 * t3 / 2;
		//position
		const double x1 = u1 + f2 * (t2 * u3 - t3 * u2) + f3 * (t1 * t2 * u2 + t1 * t3 * u3 - (t2 * t2 + t3 * t3) * u1);
		const double x2 = u2 + f2 * (t3 * u1 - t1 * u3) + f3 * (t1 * t2 * u1 + t2 * t3 * u3 - (t1 * t1 + t3 * t3) * u2);
		const double x3 = u3 + f2 * (t1 * u2 - t2 * u1) + f3 * (t1 * t3 * u1 + t2 * t3 * u2 - (t1 * t1 + t2 * t2) * u3);
		//return
		return GSE3(x1, x2, x3, q0, q1, q2, q3);
	}
	ASE3 ASE3::bracket(const ASE3& v) const
	{
		//data
		const double u1 = m_data[0];
		const double u2 = m_data[1];
		const double u3 = m_data[2];
		const double t1 = m_data[3];
		const double t2 = m_data[4];
		const double t3 = m_data[5];
		const double v1 = v.m_data[0];
		const double v2 = v.m_data[1];
		const double v3 = v.m_data[2];
		const double w1 = v.m_data[3];
		const double w2 = v.m_data[4];
		const double w3 = v.m_data[5];
		//vector
		const double r1 = t2 * w3 - t3 * w2;
		const double r2 = t3 * w1 - t1 * w3;
		const double r3 = t1 * w2 - t2 * w1;
		const double x1 = t2 * v3 - t3 * v2 - w2 * u3 + w3 * u2;
		const double x2 = t3 * v1 - t1 * v3 - w3 * u1 + w1 * u3;
		const double x3 = t1 * v2 - t2 * v1 - w1 * u2 + w2 * u1;
		//return
		return ASE3(x1, x2, x3, r1, r2, r3);
	}

	//random
	ASE3& ASE3::randu(double a, double b)
	{
		vector().randu(a, b);
		return *this;
	}

	//linear
	Vector ASE3::vector(void)
	{
		return Vector(m_data, 6);
	}
	Vector ASE3::vector(void) const
	{
		return Vector(m_data, 6);
	}
	Matrix ASE3::matrix(void) const
	{
		const double u1 = m_data[0];
		const double u2 = m_data[1];
		const double u3 = m_data[2];
		const double t1 = m_data[3];
		const double t2 = m_data[4];
		const double t3 = m_data[5];
		return Matrix({{0, -t3, t2, u1}, {t3, 0, -t1, u2}, {-t2, t1, 0, u3}, {0, 0, 0, 0}});
	}
	Matrix ASE3::adjoint(void) const
	{
		const double u1 = m_data[0];
		const double u2 = m_data[1];
		const double u3 = m_data[2];
		const double t1 = m_data[3];
		const double t2 = m_data[4];
		const double t3 = m_data[5];
		return Matrix({
			{0, -t3, t2, 0, -u3, u2},
			{t3, 0, -t1, u3, 0, -u1},
			{-t2, t1, 0, -u2, u1, 0},
			{0, 0, 0, 0, -t3, t2}, {0, 0, 0, t3, 0, -t1}, {0, 0, 0, -t2, t1, 0}
		});
	}
	Matrix ASE3::jacobian(void) const
	{
		//data
		Matrix J(6, 6);
		const Vec3 u(m_data + 0), t(m_data + 3);
		const Mat3 Su = u.spin(), St = t.spin();
		const Vec3 a = t.cross(u), b = t.cross(t.cross(u));
		//jacobian
		const double tn = t.norm();
		const double f1 = fn(tn, 1);
		const double f2 = fn(tn, 2);
		const double f3 = fn(tn, 3);
		const double f4 = fn(tn, 4);
		const double f5 = fn(tn, 5);
		const Mat3 Ri = Mat3::eye() - f1 * St + f2 * St * St;
		const Mat3 A1 = f2 * Su + f3 * (2 * St * Su - Su * St);
		const Mat3 A2 = (2 * f4 - f3) * a.outer(t) + (3 * f5 - f4) * b.outer(t);
		//jacobian
		J.span(3, 0) = 0;
		J.span(0, 3) = Ri * (A2 - A1);
		J.span(0, 0) = J.span(3, 3) = Mat3::eye() - f2 * St + f3 * St * St;
		//return
		return J;
	}
	Matrix ASE3::jacobian_inverse(void) const
	{
		//data
		Matrix Ji(6, 6);
		const Vec3 u(m_data + 0), t(m_data + 3);
		const Mat3 Su = u.spin(), St = t.spin();
		const Vec3 a = t.cross(u), b = t.cross(t.cross(u));
		//jacobian
		const double tn = t.norm();
		const double f1 = fn(tn, 1);
		const double f2 = fn(tn, 2);
		const double f3 = fn(tn, 3);
		const double f4 = fn(tn, 4);
		const double f5 = fn(tn, 5);
		const Mat3 Ri = Mat3::eye() - f1 * St + f2 * St * St;
		const Mat3 A1 = f2 * Su + f3 * (2 * St * Su - Su * St);
		const Mat3 A2 = (2 * f4 - f3) * a.outer(t) + (3 * f5 - f4) * b.outer(t);
		const Mat3 Ti = Mat3::eye() + St / 2 + (f3 - 2 * f4) / f2 / 2 * St * St;
		//jacobian
		Ji.span(3, 0) = 0;
		Ji.span(0, 0) = Ji.span(3, 3) = Ti;
		Ji.span(0, 3) = -Ti * Ri * (A2 - A1) * Ti;
		//return
		return Ji;
	}

	//operators
	ASE3& ASE3::operator*=(double s)
	{
		for(uint32_t i = 0; i < 6; i++) m_data[i] *= s;
		return *this;
	}
	ASE3& ASE3::operator/=(double s)
	{
		for(uint32_t i = 0; i < 6; i++) m_data[i] /= s;
		return *this;
	}
	ASE3& ASE3::operator+=(const ASE3& vector)
	{
		for(uint32_t i = 0; i < 6; i++) m_data[i] += vector.m_data[i];
		return *this;
	}
	ASE3& ASE3::operator-=(const ASE3& vector)
	{
		for(uint32_t i = 0; i < 6; i++) m_data[i] -= vector.m_data[i];
		return *this;
	}

	ASE3 ASE3::operator+(void) const
	{
		return *this;
	}
	ASE3 ASE3::operator-(void) const
	{
		return ASE3(*this) *= -1;
	}
	ASE3 ASE3::operator*(double s) const
	{
		return ASE3(*this) *= s;
	}
	ASE3 ASE3::operator/(double s) const
	{
		return ASE3(*this) /= s;
	}
	ASE3 ASE3::operator+(const ASE3& vector) const
	{
		return ASE3(*this) += vector;
	}
	ASE3 ASE3::operator-(const ASE3& vector) const
	{
		return ASE3(*this) -= vector;
	}

	ASE3 operator*(double s, const ASE3& vector)
	{
		return ASE3(vector) *= s;
	}

	double& ASE3::operator[](uint32_t index)
	{
		return m_data[index];
	}
	const double& ASE3::operator[](uint32_t index) const
	{
		return m_data[index];
	}

	//tests
	bool ASE3::test_adjoint(uint32_t nt, double tolerance)
	{
		//test
		ASE3 v1, v2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v1.randu(-1, +1);
			v2.randu(-1, +1);
			test = test && (v1.adjoint() * v2.vector() - v1.bracket(v2).vector()).norm() < tolerance;
		}
		//print
		printf("ASE3 adjoint test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASE3::test_bracket(uint32_t nt, double tolerance)
	{
		//test
		ASE3 v1, v2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v1.randu(-1, +1);
			v2.randu(-1, +1);
			const Matrix A1 = v1.matrix(), A2 = v2.matrix();
			test = test && (v1.bracket(v2).matrix() - (A1 * A2 - A2 * A1)).norm() < tolerance;
		}
		//print
		printf("ASE3 bracket test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASE3::test_exp_log(uint32_t nt, double tolerance)
	{
		//test
		ASE3 v;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			test = test && (v - v.exp().log()).vector().norm() < tolerance;
		}
		//print
		printf("ASE3 exp-log test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASE3::test_jacobian(uint32_t nt, double tolerance)
	{
		//test
		Matrix Jn(6, 6);
		bool test = true;
		ASE3 v, vm, vp, dm, dp;
		const double h = 1.00e-6;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			const GSE3 gi = v.exp().inverse();
			for(uint32_t j = 0; j < 6; j++)
			{
				vm = v;
				vp = v;
				vm[j] -= h;
				vp[j] += h;
				dm = (gi * vm.exp()).log();
				dp = (gi * vp.exp()).log();
				for(uint32_t k = 0; k < 6; k++) Jn(k, j) = (dp[k] - dm[k]) / (2 * h);
			}
			test = test && (Jn - v.jacobian()).norm() < tolerance;
		}
		//print
		printf("ASE3 jacobian test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASE3::test_jacobian_inverse(uint32_t nt, double tolerance)
	{
		//test
		ASE3 v;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			test = test && (v.jacobian() * v.jacobian_inverse() - Matrix::eye(6, 6)).norm() < tolerance;
		}
		//print
		printf("ASE3 jacobian inverse test:\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
}
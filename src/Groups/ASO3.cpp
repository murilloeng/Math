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
	ASO3::ASO3(const Vec3& data) : m_data{data[0], data[1], data[2]}
	{
		return;
	}
	ASO3::ASO3(const double* data) : m_data{data[0], data[1], data[2]}
	{
		return;
	}
	ASO3::ASO3(double t1, double t2, double t3) : m_data{t1, t2, t3}
	{
		return;
	}

	//group
	GSO3 ASO3::exp(void) const
	{
		const Vec3 t(m_data);
		const double tn = t.norm();
		const double f0 = fn(tn / 2, 0);
		const double f1 = fn(tn / 2, 1);
		return GSO3(f0, f1 * t[0] / 2, f1 * t[1] / 2, f1 * t[2] / 2);
	}
	ASO3 ASO3::bracket(const ASO3& v) const
	{
		return ASO3(Vec3(m_data).cross(Vec3(v.m_data)));
	}

	//random
	ASO3& ASO3::randu(double a, double b)
	{
		return *this = ASO3(Vec3().randu(a, b).data());
	}

	//linear
	Vec3 ASO3::vector(void) const
	{
		return Vec3(m_data);
	}
	Mat3 ASO3::matrix(void) const
	{
		return Vec3(m_data).spin();
	}
	Mat3 ASO3::adjoint(void) const
	{
		return Vec3(m_data).spin();
	}
	Mat3 ASO3::jacobian(void) const
	{
		const Vec3 t(m_data);
		const Mat3 T = t.spin();
		const double tn = t.norm();
		const double f2 = fn(tn, 2);
		const double f3 = fn(tn, 3);
		return Mat3(Mode::eye) - f2 * T + f3 * T * T;
	}
	Mat3 ASO3::jacobian_inverse(void) const
	{
		const Vec3 t(m_data);
		const Mat3 T = t.spin();
		const double tn = t.norm();
		const double f2 = fn(tn, 2);
		const double f3 = fn(tn, 3);
		const double f4 = fn(tn, 4);
		return Mat3(Mode::eye) + T / 2 + (f3 - 2 * f4) / f2 / 2 * T * T;
	}
	Mat3 ASO3::jacobian_gradient(const Vec3& vector) const
	{
		//data
		const Vec3 t(m_data);
		const Vec3 u(vector);
		const double tn = t.norm();
		const double f2 = fn(tn, 2);
		const double f3 = fn(tn, 3);
		const double f4 = fn(tn, 4);
		const double f5 = fn(tn, 5);
		//gradient
		Mat3 G(Mode::zeros);
		G += (2 * f4 - f3) * t.cross(u).outer(t);
		G += (3 * f5 - f4) * t.cross(t.cross(u)).outer(t);
		G -= f2 * u.spin() + f3 * t.cross(u).spin() + f3 * t.spin() * u.spin();
		//return
		return G;
	}
	Mat3 ASO3::jacobian_inverse_gradient(const Vec3& v, bool transpose) const
	{
		const Vec3 t(m_data);
		const Mat3 T = t.spin();
		const Mat3 V = v.spin();
		const Mat3 M = v.outer(t);
		const double tn = t.norm();
		const double f2 = fn(tn, 2);
		const double f3 = fn(tn, 3);
		const double f4 = fn(tn, 4);
		const double f5 = fn(tn, 5);
		const double f6 = fn(tn, 6);
		const double st = transpose ? +1 : -1;
		return st * V / 2 + (f3 - 2 * f4) / f2 / 2 * (V * T - 2 * T * V) + (f5 - 4 * f6) / f2 / 2 * T * T * M;
	}

	//operators
	ASO3& ASO3::operator*=(double s)
	{
		for(uint32_t i = 0; i < 3; i++) m_data[i] *= s;
		return *this;
	}
	ASO3& ASO3::operator/=(double s)
	{
		for(uint32_t i = 0; i < 3; i++) m_data[i] /= s;
		return *this;
	}
	ASO3& ASO3::operator+=(const ASO3& vector)
	{
		for(uint32_t i = 0; i < 3; i++) m_data[i] += vector.m_data[i];
		return *this;
	}
	ASO3& ASO3::operator-=(const ASO3& vector)
	{
		for(uint32_t i = 0; i < 3; i++) m_data[i] -= vector.m_data[i];
		return *this;
	}

	ASO3 ASO3::operator+(void) const
	{
		return *this;
	}
	ASO3 ASO3::operator-(void) const
	{
		return ASO3(*this) *= -1;
	}
	ASO3 ASO3::operator*(double s) const
	{
		return ASO3(*this) *= s;
	}
	ASO3 ASO3::operator/(double s) const
	{
		return ASO3(*this) /= s;
	}
	ASO3 ASO3::operator+(const ASO3& vector) const
	{
		return ASO3(*this) += vector;
	}
	ASO3 ASO3::operator-(const ASO3& vector) const
	{
		return ASO3(*this) -= vector;
	}

	ASO3 operator*(double s, const ASO3& vector)
	{
		return ASO3(vector) *= s;
	}

	double& ASO3::operator[](uint32_t index)
	{
		return m_data[index];
	}
	const double& ASO3::operator[](uint32_t index) const
	{
		return m_data[index];
	}

	//tests
	bool ASO3::test_adjoint(uint32_t nt, double tolerance)
	{
		//test
		ASO3 v1, v2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v1.randu(-1, +1);
			v2.randu(-1, +1);
			test = test && (v1.adjoint() * v2.vector() - v1.bracket(v2).vector()).norm() < tolerance;
		}
		//print
		printf("ASO3 adjoint test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASO3::test_bracket(uint32_t nt, double tolerance)
	{
		//test
		ASO3 v1, v2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v1.randu(-1, +1);
			v2.randu(-1, +1);
			const Mat3 A1 = v1.matrix(), A2 = v2.matrix();
			test = test && (v1.bracket(v2).matrix() - (A1 * A2 - A2 * A1)).norm() < tolerance;
		}
		//print
		printf("ASO3 bracket test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASO3::test_exp_log(uint32_t nt, double tolerance)
	{
		//test
		ASO3 v;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			test = test && (v - v.exp().log()).vector().norm() < tolerance;
		}
		//print
		printf("ASO3 exp-log test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASO3::test_jacobian(uint32_t nt, double tolerance)
	{
		//test
		Mat3 Jn;
		bool test = true;
		ASO3 v, vm, vp, dm, dp;
		const double h = 1.00e-6;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			const GSO3 gi = v.exp().inverse();
			for(uint32_t j = 0; j < 3; j++)
			{
				vm = v;
				vp = v;
				vm[j] -= h;
				vp[j] += h;
				dm = (gi * vm.exp()).log();
				dp = (gi * vp.exp()).log();
				for(uint32_t k = 0; k < 3; k++) Jn(k, j) = (dp[k] - dm[k]) / (2 * h);
			}
			test = test && (Jn - v.jacobian()).norm() < tolerance;
		}
		//print
		printf("ASO3 jacobian test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASO3::test_jacobian_inverse(uint32_t nt, double tolerance)
	{
		//test
		ASO3 v;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			test = test && (v.jacobian() * v.jacobian_inverse() - Mat3(Mode::eye)).norm() < tolerance;
		}
		//print
		printf("ASO3 jacobian inverse test:\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
}
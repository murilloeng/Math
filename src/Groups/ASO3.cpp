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
		return GSO3(m_data[0], m_data[1], m_data[2]);
	}
	ASO3 ASO3::bracket(const ASO3& v) const
	{
		const double t1 = m_data[0];
		const double t2 = m_data[1];
		const double t3 = m_data[2];
		const double v1 = v.m_data[0];
		const double v2 = v.m_data[1];
		const double v3 = v.m_data[2];
		return ASO3(t2 * v3 - t3 * v2, t3 * v1 - t1 * v3, t1 * v2 - t2 * v1);
	}

	//random
	ASO3& ASO3::randu(double a, double b)
	{
		vector().randu(a, b);
		return *this;
	}

	//linear
	Vector ASO3::vector(void)
	{
		return Vector(m_data, 3);
	}
	Vector ASO3::vector(void) const
	{
		return Vector(m_data, 3);
	}
	Matrix ASO3::matrix(void) const
	{
		const double t1 = m_data[0];
		const double t2 = m_data[1];
		const double t3 = m_data[2];
		return Matrix({{0, -t3, t2}, {t3, 0, -t1}, {-t2, t1, 0}});
	}
	Matrix ASO3::adjoint(void) const
	{
		const double t1 = m_data[0];
		const double t2 = m_data[1];
		const double t3 = m_data[2];
		return Matrix({{0, -t3, t2}, {t3, 0, -t1}, {-t2, t1, 0}});
	}
	Matrix ASO3::jacobian(void) const
	{
		const Matrix A = matrix();
		const double t1 = m_data[0];
		const double t2 = m_data[1];
		const double t3 = m_data[2];
		const double tn = sqrt(t1 * t1 + t2 * t2 + t3 * t3);
		return Matrix::eye(3, 3) - fn(tn, 2) * A + fn(tn, 3) * A * A;
	}
	Matrix ASO3::jacobian_inverse(void) const
	{
		const Matrix A = matrix();
		const double t1 = m_data[0];
		const double t2 = m_data[1];
		const double t3 = m_data[2];
		const double tn = sqrt(t1 * t1 + t2 * t2 + t3 * t3);
		return Matrix::eye(3, 3) + A / 2 + (fn(tn, 3) - 2 * fn(tn, 4)) / fn(tn, 2) / 2 * A * A;
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
			const Matrix A1 = v1.matrix(), A2 = v2.matrix();
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
		Matrix Jn(3, 3);
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
			test = test && (v.jacobian() * v.jacobian_inverse() - Matrix::eye(3, 3)).norm() < tolerance;
		}
		//print
		printf("ASO3 jacobian inverse test:\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
}
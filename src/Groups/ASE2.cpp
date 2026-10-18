//std
#include <cstdio>

//Math
#include "Math/inc/Groups/ASE2.hpp"
#include "Math/inc/Groups/GSE2.hpp"
#include "Math/inc/Miscellaneous/rotations.hpp"

namespace math
{
	//constructor
	ASE2::ASE2(const double* data) : m_data{data[0], data[1], data[2]}
	{
		return;
	}
	ASE2::ASE2(double x1, double x2, double t) : m_data{x1, x2, t}
	{
		return;
	}

	//group
	GSE2 ASE2::exp(void) const
	{
		const double u1 = m_data[0];
		const double u2 = m_data[1];
		const double tr = m_data[2];
		const double sr = fn(tr, 1);
		const double cr = fn(tr, 2) * tr;
		return GSE2(sr * u1 - cr * u2, sr * u2 + cr * u1, tr);
	}
	ASE2 ASE2::bracket(const ASE2& vector) const
	{
		const double u1 = m_data[0];
		const double u2 = m_data[1];
		const double t1 = m_data[2];
		const double v1 = vector.m_data[0];
		const double v2 = vector.m_data[1];
		const double t2 = vector.m_data[2];
		return ASE2(t2 * u2 - t1 * v2, t1 * v1 - t2 * u1, 0);
	}

	//random
	ASE2& ASE2::randu(double a, double b)
	{
		vector().randu(a, b);
		return *this;
	}

	//linear
	Vector ASE2::vector(void)
	{
		return Vector(m_data, 3);
	}
	Vector ASE2::vector(void) const
	{
		return Vector(m_data, 3);
	}
	Matrix ASE2::matrix(void) const
	{
		const double u1 = m_data[0];
		const double u2 = m_data[1];
		const double tr = m_data[2];
		return Matrix({{0, -tr, u1}, {tr, 0, u2}, {0, 0, 0}});
	}
	Matrix ASE2::adjoint(void) const
	{
		const double u1 = m_data[0];
		const double u2 = m_data[1];
		const double tr = m_data[2];
		return Matrix({{0, -tr, u2}, {tr, 0, -u1}, {0, 0, 0}});
	}
	Matrix ASE2::jacobian(void) const
	{
		const double u1 = m_data[0];
		const double u2 = m_data[1];
		const double tr = m_data[2];
		const double f1 = fn(tr, 1);
		const double f2 = fn(tr, 2);
		const double f3 = fn(tr, 3);
		return Matrix({
			{f1, +tr * f2, tr * f3 * u1 - f2 * u2},
			{-tr * f2, f1, tr * f3 * u2 + f2 * u1},
			{0, 0, 1}
		});
	}
	Matrix ASE2::jacobian_inverse(void) const
	{
		const double u1 = m_data[0];
		const double u2 = m_data[1];
		const double tr = m_data[2];
		const double f1 = fn(tr, 1);
		const double f2 = fn(tr, 2);
		const double f3 = fn(tr, 3);
		const double f4 = fn(tr, 4);
		return Matrix({
			{f1 / f2 / 2, -tr / 2, tr * (f3 - 2 * f4) / f2 / 2 * u1 + u2 / 2},
			{+tr / 2, f1 / f2 / 2, tr * (f3 - 2 * f4) / f2 / 2 * u2 - u1 / 2},
			{0, 0, 1}
		});
	}

	//operators
	ASE2& ASE2::operator*=(double s)
	{
		for(uint32_t i = 0; i < 3; i++) m_data[i] *= s;
		return *this;
	}
	ASE2& ASE2::operator/=(double s)
	{
		for(uint32_t i = 0; i < 3; i++) m_data[i] /= s;
		return *this;
	}
	ASE2& ASE2::operator+=(const ASE2& vector)
	{
		for(uint32_t i = 0; i < 3; i++) m_data[i] += vector.m_data[i];
		return *this;
	}
	ASE2& ASE2::operator-=(const ASE2& vector)
	{
		for(uint32_t i = 0; i < 3; i++) m_data[i] -= vector.m_data[i];
		return *this;
	}

	ASE2 ASE2::operator+(void) const
	{
		return *this;
	}
	ASE2 ASE2::operator-(void) const
	{
		return ASE2(*this) *= -1;
	}
	ASE2 ASE2::operator*(double s) const
	{
		return ASE2(*this) *= s;
	}
	ASE2 ASE2::operator/(double s) const
	{
		return ASE2(*this) /= s;
	}
	ASE2 ASE2::operator+(const ASE2& vector) const
	{
		return ASE2(*this) += vector;
	}
	ASE2 ASE2::operator-(const ASE2& vector) const
	{
		return ASE2(*this) -= vector;
	}

	ASE2 operator*(double s, const ASE2& vector)
	{
		return ASE2(vector) *= s;
	}

	double& ASE2::operator[](uint32_t index)
	{
		return m_data[index];
	}
	const double& ASE2::operator[](uint32_t index) const
	{
		return m_data[index];
	}

	//tests
	bool ASE2::test_adjoint(uint32_t nt, double tolerance)
	{
		//test
		ASE2 v1, v2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v1.vector().randu(-1, +1);
			v2.vector().randu(-1, +1);
			test = test && (v1.adjoint() * v2.vector() - v1.bracket(v2).vector()).norm() < tolerance;
		}
		//print
		printf("ASE2 adjoint test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASE2::test_bracket(uint32_t nt, double tolerance)
	{
		//test
		ASE2 v1, v2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v1.vector().randu();
			v2.vector().randu();
			const Matrix A1 = v1.matrix(), A2 = v2.matrix();
			test = test && (v1.bracket(v2).matrix() - (A1 * A2 - A2 * A1)).norm() < tolerance;
		}
		//print
		printf("ASE2 bracket test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASE2::test_exp_log(uint32_t nt, double tolerance)
	{
		//test
		ASE2 v;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			test = test && (v - v.exp().log()).vector().norm() < tolerance;
		}
		//print
		printf("ASE2 exp-log test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASE2::test_jacobian(uint32_t nt, double tolerance)
	{
		//test
		Matrix Jn(3, 3);
		bool test = true;
		ASE2 v, vm, vp, dm, dp;
		const double h = 1.00e-6;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			const GSE2 gi = v.exp().inverse();
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
		printf("ASE2 jacobian test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASE2::test_jacobian_inverse(uint32_t nt, double tolerance)
	{
		//test
		ASE2 v;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			test = test && (v.jacobian() * v.jacobian_inverse() - Matrix::eye(3, 3)).norm() < tolerance;
		}
		//print
		printf("ASE2 jacobian inverse test:\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
}
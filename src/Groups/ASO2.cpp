//Math
#include "Math/inc/Groups/ASO2.hpp"
#include "Math/inc/Groups/GSO2.hpp"
#include "Math/inc/Miscellaneous/util.hpp"

namespace math
{
	//constructor
	ASO2::ASO2(double data) : m_data{data}
	{
		return;
	}
	ASO2::ASO2(const double* data) : m_data{data[0]}
	{
		return;
	}

	//group
	GSO2 ASO2::exp(void) const
	{
		return GSO2(m_data);
	}
	ASO2 ASO2::bracket(const ASO2&) const
	{
		return ASO2(0.0);
	}

	//random
	ASO2& ASO2::randu(double a, double b)
	{
		vector().randu(a, b);
		return *this;
	}

	//linear
	Vector ASO2::vector(void)
	{
		return Vector(&m_data, 1);
	}
	Vector ASO2::vector(void) const
	{
		return Vector(&m_data, 1);
	}
	Matrix ASO2::matrix(void) const
	{
		return Matrix({{0, -m_data}, {m_data, 0}});
	}
	Matrix ASO2::adjoint(void) const
	{
		return Matrix({{0}});
	}
	Matrix ASO2::jacobian(void) const
	{
		return Matrix({{1}});
	}
	Matrix ASO2::jacobian_inverse(void) const
	{
		return Matrix({{1}});
	}
	Matrix ASO2::jacobian_gradient(const Vector&) const
	{
		return Matrix({{0}});
	}

	//operators
	ASO2& ASO2::operator*=(double s)
	{
		m_data *= s;
		return *this;
	}
	ASO2& ASO2::operator/=(double s)
	{
		m_data /= s;
		return *this;
	}
	ASO2& ASO2::operator+=(const ASO2& vector)
	{
		m_data += vector.m_data;
		return *this;
	}
	ASO2& ASO2::operator-=(const ASO2& vector)
	{
		m_data -= vector.m_data;
		return *this;
	}

	ASO2 ASO2::operator+(void) const
	{
		return *this;
	}
	ASO2 ASO2::operator-(void) const
	{
		return ASO2(*this) *= -1;
	}
	ASO2 ASO2::operator*(double s) const
	{
		return ASO2(*this) *= s;
	}
	ASO2 ASO2::operator/(double s) const
	{
		return ASO2(*this) /= s;
	}
	ASO2 ASO2::operator+(const ASO2& vector) const
	{
		return ASO2(*this) += vector;
	}
	ASO2 ASO2::operator-(const ASO2& vector) const
	{
		return ASO2(*this) -= vector;
	}

	ASO2 operator*(double s, const ASO2& vector)
	{
		return ASO2(s * vector.m_data);
	}

	double& ASO2::operator[](uint32_t)
	{
		return m_data;
	}
	const double& ASO2::operator[](uint32_t) const
	{
		return m_data;
	}

	//tests
	bool ASO2::test_adjoint(uint32_t nt, double tolerance)
	{
		//test
		ASO2 v1, v2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v1.randu(-1, +1);
			v2.randu(-1, +1);
			test = test && (v1.adjoint() * v2.vector() - v1.bracket(v2).vector()).norm() < tolerance;
		}
		//print
		printf("ASO2 adjoint test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASO2::test_bracket(uint32_t nt, double tolerance)
	{
		//test
		ASO2 v1, v2;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v1.randu(-1, +1);
			v2.randu(-1, +1);
			const Matrix A1 = v1.matrix(), A2 = v2.matrix();
			test = test && (v1.bracket(v2).matrix() - (A1 * A2 - A2 * A1)).norm() < tolerance;
		}
		//print
		printf("ASO2 bracket test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASO2::test_exp_log(uint32_t nt, double tolerance)
	{
		//test
		ASO2 v;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			test = test && (v - v.exp().log()).vector().norm() < tolerance;
		}
		//print
		printf("ASO2 exp-log test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASO2::test_jacobian(uint32_t nt, double tolerance)
	{
		//test
		ASO2 v;
		Matrix Jn(1, 1);
		bool test = true;
		const double h = 1.00e-06;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			const GSO2 gi = v.exp().inverse();
			const ASO2 vm(v.m_data - h), vp(v.m_data + h);
			Jn[0] = ((gi * vp.exp()).log().m_data - (gi * vm.exp()).log().m_data) / (2 * h);
			test = test && (Jn - v.jacobian()).norm() < tolerance;
		}
		//print
		printf("ASO2 jacobian test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool ASO2::test_jacobian_inverse(uint32_t nt, double tolerance)
	{
		//test
		ASO2 v;
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			v.randu(-1, +1);
			test = test && (v.jacobian() * v.jacobian_inverse() - Matrix::eye(1, 1)).norm() < tolerance;
		}
		//print
		printf("ASO2 jacobian inverse test:\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
}
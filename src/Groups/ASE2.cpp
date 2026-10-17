//Math
#include "Math/inc/Groups/ASE2.hpp"
#include "Math/inc/Groups/GSE2.hpp"

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
		return GSE2(m_data);
	}
	ASE2 ASE2::bracket(const ASE2& vector) const
	{
		const double x1 = m_data[0];
		const double y1 = m_data[1];
		const double t1 = m_data[2];
		const double x2 = vector.m_data[0];
		const double y2 = vector.m_data[1];
		const double t2 = vector.m_data[2];
		return ASE2(t2 * y1 - t1 * y2, t1 * x2 - t2 * x1, 0);
	}

	//linear
	Vector ASE2::vector(void) const
	{
		return Vector(m_data, 3);
	}
	Matrix ASE2::matrix(void) const
	{
		const double x = m_data[0];
		const double y = m_data[1];
		const double t = m_data[2];
		return Matrix({{0, -t, x}, {t, 0, y}, {0, 0, 0}});
	}
	Matrix ASE2::adjoint(void) const
	{
		const double x = m_data[0];
		const double y = m_data[1];
		const double t = m_data[2];
		return Matrix({{0, -t, y}, {t, 0, -x}, {0, 0, 0}});
	}
	Matrix ASE2::gradient(void) const
	{
		return Matrix({{1}});
	}
	Matrix ASE2::gradient_inverse(void) const
	{
		return Matrix({{1}});
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
}
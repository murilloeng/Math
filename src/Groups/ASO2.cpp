//Math
#include "Math/inc/Groups/ASO2.hpp"
#include "Math/inc/Groups/GSO2.hpp"

namespace math
{
	//constructor
	ASO2::ASO2(double data) : m_data{data}
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
		return ASO2(0);
	}

	//linear
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
}
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

	//linear
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
}
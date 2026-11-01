//std
#include <cmath>

//Math
#include "Math/inc/Linear2/Vec2.hpp"
#include "Math/inc/Linear2/Mat2.hpp"

namespace math
{
	namespace linear
	{
		//constructors
		Vec2::Vec2(const double* data) : m_data{data[0], data[1]}
		{
			return;
		}
		Vec2::Vec2(double v0, double v1) : m_data{v0, v1}
		{
			return;
		}

		//data
		double* Vec2::data(void)
		{
			return m_data;
		}
		const double* Vec2::data(void) const
		{
			return m_data;
		}

		//operators
		Vec2 Vec2::operator+(void) const
		{
			return *this;
		}
		Vec2 Vec2::operator-(void) const
		{
			return Vec2(*this) *= -1;
		}

		Vec2 Vec2::operator/(double s) const
		{
			return Vec2(*this) /= s;
		}
		Vec2 Vec2::operator+(const Vec2& v) const
		{
			return Vec2(*this) += v;
		}
		Vec2 Vec2::operator-(const Vec2& v) const
		{
			return Vec2(*this) -= v;
		}

		Vec2& Vec2::operator+=(double s)
		{
			m_data[0] += s;
			m_data[1] += s;
			return *this;
		}
		Vec2& Vec2::operator-=(double s)
		{
			m_data[0] -= s;
			m_data[1] -= s;
			return *this;
		}
		Vec2& Vec2::operator*=(double s)
		{
			m_data[0] *= s;
			m_data[1] *= s;
			return *this;
		}
		Vec2& Vec2::operator/=(double s)
		{
			m_data[0] /= s;
			m_data[1] /= s;
			return *this;
		}

		Vec2& Vec2::operator+=(const Vec2& v)
		{
			m_data[0] += v.m_data[0];
			m_data[1] += v.m_data[1];
			return *this;
		}
		Vec2& Vec2::operator-=(const Vec2& v)
		{
			m_data[0] -= v.m_data[0];
			m_data[1] -= v.m_data[1];
			return *this;
		}

		double& Vec2::operator[](uint32_t i)
		{
			return m_data[i];
		}
		const double& Vec2::operator[](uint32_t i) const
		{
			return m_data[i];
		}

		//norm
		double Vec2::norm(void) const
		{
			return sqrt(inner(*this));
		}

		//direction
		Vec2 Vec2::unit(void) const
		{
			return *this / norm();
		}
		Vec2& Vec2::normalize(void)
		{
			return *this /= norm();
		}

		//projection
		Mat2 Vec2::projection(void) const
		{
			const double v0 = m_data[0];
			const double v1 = m_data[1];
			return {{1 - v0 * v0, -v0 * v1}, {-v0 * v1, 1 - v1 * v1}};
		}

		//products
		double Vec2::inner(const Vec2& v) const
		{
			return m_data[0] * v.m_data[0] + m_data[1] * v.m_data[1];
		}
		double Vec2::cross(const Vec2& v) const
		{
			return m_data[0] * v.m_data[1] - m_data[1] * v.m_data[0];
		}

		//interpolation
		Vec2 Vec2::lerp(const Vec2& v, double s) const
		{
			Vec2 r;
			r.m_data[0] = (1 - s) * m_data[0] + s * v.m_data[0];
			r.m_data[1] = (1 - s) * m_data[1] + s * v.m_data[1];
			return r;
		}

		//friends
		Vec2 operator*(double s, const Vec2& v)
		{
			return Vec2(v) *= s;
		}
	}
}
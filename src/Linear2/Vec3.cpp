//std
#include <cmath>

//Math
#include "Math/inc/Linear2/Vec3.hpp"
#include "Math/inc/Linear2/Mat3.hpp"

namespace math
{
	namespace linear
	{
		//constructors
		Vec3::Vec3(const double* data) : m_data{data[0], data[1], data[2]}
		{
			return;
		}
		Vec3::Vec3(double v0, double v1, double v2) : m_data{v0, v1, v2}
		{
			return;
		}

		//data
		double* Vec3::data(void)
		{
			return m_data;
		}
		const double* Vec3::data(void) const
		{
			return m_data;
		}

		//operators
		Vec3 Vec3::operator+(void) const
		{
			return *this;
		}
		Vec3 Vec3::operator-(void) const
		{
			return Vec3(*this) *= -1;
		}

		Vec3 Vec3::operator/(double s) const
		{
			return Vec3(*this) /= s;
		}
		Vec3 Vec3::operator+(const Vec3& v) const
		{
			return Vec3(*this) += v;
		}
		Vec3 Vec3::operator-(const Vec3& v) const
		{
			return Vec3(*this) -= v;
		}

		Vec3& Vec3::operator+=(double s)
		{
			m_data[0] += s;
			m_data[1] += s;
			m_data[2] += s;
			return *this;
		}
		Vec3& Vec3::operator-=(double s)
		{
			m_data[0] -= s;
			m_data[1] -= s;
			m_data[2] -= s;
			return *this;
		}
		Vec3& Vec3::operator*=(double s)
		{
			m_data[0] *= s;
			m_data[1] *= s;
			m_data[2] *= s;
			return *this;
		}
		Vec3& Vec3::operator/=(double s)
		{
			m_data[0] /= s;
			m_data[1] /= s;
			m_data[2] /= s;
			return *this;
		}

		Vec3& Vec3::operator+=(const Vec3& v)
		{
			m_data[0] += v.m_data[0];
			m_data[1] += v.m_data[1];
			m_data[2] += v.m_data[2];
			return *this;
		}
		Vec3& Vec3::operator-=(const Vec3& v)
		{
			m_data[0] -= v.m_data[0];
			m_data[1] -= v.m_data[1];
			m_data[2] -= v.m_data[2];
			return *this;
		}

		double& Vec3::operator[](uint32_t i)
		{
			return m_data[i];
		}
		const double& Vec3::operator[](uint32_t i) const
		{
			return m_data[i];
		}

		//norm
		double Vec3::norm(void) const
		{
			return sqrt(inner(*this));
		}

		//direction
		Vec3 Vec3::unit(void) const
		{
			return *this / norm();
		}
		Vec3& Vec3::normalize(void)
		{
			return *this /= norm();
		}
		Vec3& Vec3::project(const Vec3& v)
		{
			return *this -= inner(v) * v;
		}

		//projection
		Mat3 Vec3::spin(void) const
		{
			Mat3 r;
			r.m_data[7] = -m_data[0];
			r.m_data[2] = -m_data[1];
			r.m_data[3] = -m_data[2];
			r.m_data[5] = +m_data[0];
			r.m_data[6] = +m_data[1];
			r.m_data[1] = +m_data[2];
			r.m_data[0] = r.m_data[4] = r.m_data[8] = 0;
			return r;
		}
		Mat3 Vec3::projection(void) const
		{
			Mat3 r;
			r.m_data[0] = 1 - m_data[0] * m_data[0];
			r.m_data[4] = 1 - m_data[1] * m_data[1];
			r.m_data[8] = 1 - m_data[2] * m_data[2];
			r.m_data[1] = r.m_data[3] = -m_data[0] * m_data[1];
			r.m_data[2] = r.m_data[6] = -m_data[0] * m_data[2];
			r.m_data[5] = r.m_data[7] = -m_data[1] * m_data[2];
			return r;
		}

		Mat3 Vec3::outer(const Vec3& v) const
		{
			Mat3 r;
			r.m_data[0] = m_data[0] * v.m_data[0];
			r.m_data[1] = m_data[1] * v.m_data[0];
			r.m_data[2] = m_data[2] * v.m_data[0];
			r.m_data[3] = m_data[0] * v.m_data[1];
			r.m_data[4] = m_data[1] * v.m_data[1];
			r.m_data[5] = m_data[2] * v.m_data[1];
			r.m_data[6] = m_data[0] * v.m_data[2];
			r.m_data[7] = m_data[1] * v.m_data[2];
			r.m_data[8] = m_data[2] * v.m_data[2];
			return r;
		}
		Vec3 Vec3::cross(const Vec3& v) const
		{
			Vec3 r;
			r.m_data[0] = m_data[1] * v.m_data[2] - m_data[2] * v.m_data[1];
			r.m_data[1] = m_data[2] * v.m_data[0] - m_data[0] * v.m_data[2];
			r.m_data[2] = m_data[0] * v.m_data[1] - m_data[1] * v.m_data[0];
			return r;
		}
		double Vec3::inner(const Vec3& v) const
		{
			return m_data[0] * v.m_data[0] + m_data[1] * v.m_data[1] + m_data[2] * v.m_data[2];
		}

		//interpolation
		Vec3 Vec3::lerp(const Vec3& v, double s) const
		{
			Vec3 r;
			r.m_data[0] = (1 - s) * m_data[0] + s * v.m_data[0];
			r.m_data[1] = (1 - s) * m_data[1] + s * v.m_data[1];
			r.m_data[2] = (1 - s) * m_data[2] + s * v.m_data[2];
			return r;
		}

		//friends
		Vec3 operator*(double s, const Vec3& v)
		{
			return Vec3(v) *= s;
		}
	}
}
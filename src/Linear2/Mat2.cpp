//std
#include <cmath>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear2/Vec2.hpp"
#include "Math/inc/Linear2/Mat2.hpp"

namespace math
{
	namespace linear
	{
		//constructors
		Mat2::Mat2(void)
		{
			memset(m_data, 0, 4 * sizeof(double));
		}
		Mat2::Mat2(const double* data)
		{
			memcpy(m_data, data, 4 * sizeof(double));
		}
		Mat2::Mat2(const Vec2& v1, const Vec2& v2)
		{
			memcpy(m_data + 0, v1.m_data, 2 * sizeof(double));
			memcpy(m_data + 2, v2.m_data, 2 * sizeof(double));
		}
		Mat2::Mat2(std::initializer_list<std::initializer_list<double>> data)
		{
			const std::initializer_list<double>* ptr = data.begin();
			if(data.size() != 2 || ptr[0].size() != 2 || ptr[1].size() != 2)
			{
				throw std::runtime_error("Mat2 constructor has incompatible dimensions!");
			}
			memcpy(m_data + 0, ptr[0].begin(), 2 * sizeof(double));
			memcpy(m_data + 2, ptr[1].begin(), 2 * sizeof(double));
		}

		//data
		double* Mat2::data(void)
		{
			return m_data;
		}
		const double* Mat2::data(void) const
		{
			return m_data;
		}

		//operators
		Mat2 Mat2::operator+(void) const
		{
			return *this;
		}
		Mat2 Mat2::operator-(void) const
		{
			return Mat2(*this) *= -1;
		}

		Mat2 Mat2::operator/(double s) const
		{
			return Mat2(*this) /= s;
		}
		Mat2 Mat2::operator+(const Mat2& m) const
		{
			return Mat2(*this) += m;
		}
		Mat2 Mat2::operator-(const Mat2& m) const
		{
			return Mat2(*this) -= m;
		}
		Mat2 Mat2::operator*(const Mat2& m) const
		{
			Mat2 r;
			r.m_data[0] = m_data[0] * m.m_data[0] + m_data[2] * m.m_data[1];
			r.m_data[1] = m_data[1] * m.m_data[0] + m_data[3] * m.m_data[1];
			r.m_data[2] = m_data[0] * m.m_data[2] + m_data[2] * m.m_data[3];
			r.m_data[3] = m_data[1] * m.m_data[2] + m_data[3] * m.m_data[3];
			return r;
		}
		Vec2 Mat2::operator*(const Vec2& v) const
		{
			Vec2 r;
			r[0] = m_data[0] * v.m_data[0] + m_data[2] * v.m_data[1];
			r[1] = m_data[1] * v.m_data[0] + m_data[3] * v.m_data[1];
			return r;
		}

		Mat2& Mat2::operator+=(double s)
		{
			for(uint32_t i = 0; i < 4; i++)
			{
				m_data[i] += s;
			}
			return *this;
		}
		Mat2& Mat2::operator-=(double s)
		{
			for(uint32_t i = 0; i < 4; i++)
			{
				m_data[i] -= s;
			}
			return *this;
		}
		Mat2& Mat2::operator*=(double s)
		{
			for(uint32_t i = 0; i < 4; i++)
			{
				m_data[i] *= s;
			}
			return *this;
		}
		Mat2& Mat2::operator/=(double s)
		{
			for(uint32_t i = 0; i < 4; i++)
			{
				m_data[i] /= s;
			}
			return *this;
		}

		Mat2& Mat2::operator+=(const Mat2& m)
		{
			for(uint32_t i = 0; i < 4; i++)
			{
				m_data[i] += m.m_data[i];
			}
			return *this;
		}
		Mat2& Mat2::operator-=(const Mat2& m)
		{
			for(uint32_t i = 0; i < 4; i++)
			{
				m_data[i] -= m.m_data[i];
			}
			return *this;
		}

		double& Mat2::operator[](uint32_t i)
		{
			return m_data[i];
		}
		double& Mat2::operator()(uint32_t i, uint32_t j)
		{
			return m_data[i + 2 * j];
		}

		const double& Mat2::operator[](uint32_t i) const
		{
			return m_data[i];
		}
		const double& Mat2::operator()(uint32_t i, uint32_t j) const
		{
			return m_data[i + 2 * j];
		}

		//linear
		Mat2 Mat2::eye(void)
		{
			return Mat2({{1, 0}, {0, 1}});
		}
		Mat2 Mat2::inverse(void) const
		{
			Mat2 r;
			const double d = determinant();
			r.m_data[0] = +m_data[3] / d;
			r.m_data[1] = -m_data[1] / d;
			r.m_data[2] = -m_data[2] / d;
			r.m_data[3] = +m_data[0] / d;
			return r;
		}
		Mat2 Mat2::transpose(void) const
		{
			Mat2 r(*this);
			r.m_data[1] = m_data[2];
			r.m_data[2] = m_data[1];
			return r;
		}
		Vec2 Mat2::eigenvalues(void) const
		{
			Vec2 s;
			const double a = m_data[0] + m_data[3];
			const double b = m_data[0] - m_data[3];
			const double c = m_data[1] * m_data[2];
			s[0] = a / 2 - sqrt(b * b + 4 * c) / 2;
			s[1] = a / 2 + sqrt(b * b + 4 * c) / 2;
			return s;
		}

		double Mat2::trace(void) const
		{
			return m_data[0] + m_data[3];
		}
		double Mat2::determinant(void) const
		{
			return m_data[0] * m_data[3] - m_data[1] * m_data[2];
		}
		double Mat2::invariant_I1(void) const
		{
			return m_data[0] + m_data[3];
		}
		double Mat2::invariant_I2(void) const
		{
			return m_data[0] * m_data[3] - m_data[1] * m_data[2];
		}
		double Mat2::deviatoric_J2(void) const
		{
			return (m_data[0] - m_data[3]) * (m_data[0] - m_data[3]) / 4 + m_data[1] * m_data[2];
		}

		//friends
		Mat2 operator*(double s, const Mat2& m)
		{
			return Mat2(m) *= s;
		}
	}
}
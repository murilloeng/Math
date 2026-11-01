//std
#include <cmath>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear2/Vec3.hpp"
#include "Math/inc/Linear2/Mat3.hpp"

namespace math
{
	namespace linear
	{
		//constructors
		Mat3::Mat3(void)
		{
			memset(m_data, 0, 9 * sizeof(double));
		}
		Mat3::Mat3(const double* data)
		{
			memcpy(m_data, data, 9 * sizeof(double));
		}
		Mat3::Mat3(const Vec3& v1, const Vec3& v2, const Vec3& v3)
		{
			memcpy(m_data + 0, v1.m_data, 3 * sizeof(double));
			memcpy(m_data + 3, v2.m_data, 3 * sizeof(double));
			memcpy(m_data + 6, v3.m_data, 3 * sizeof(double));
		}
		Mat3::Mat3(std::initializer_list<std::initializer_list<double>> data)
		{
			const std::initializer_list<double>* ptr = data.begin();
			if(data.size() != 3 || ptr[0].size() != 3 || ptr[1].size() != 3 || ptr[2].size() != 3)
			{
				throw std::runtime_error("Mat3 constructor has incompatible dimensions!");
			}
			memcpy(m_data + 0, ptr[0].begin(), 3 * sizeof(double));
			memcpy(m_data + 3, ptr[1].begin(), 3 * sizeof(double));
			memcpy(m_data + 6, ptr[2].begin(), 3 * sizeof(double));
		}

		//data
		double* Mat3::data(void)
		{
			return m_data;
		}
		const double* Mat3::data(void) const
		{
			return m_data;
		}

		//operators
		Mat3 Mat3::operator+(void) const
		{
			return *this;
		}
		Mat3 Mat3::operator-(void) const
		{
			return Mat3(*this) *= -1;
		}

		Mat3 Mat3::operator/(double s) const
		{
			return Mat3(*this) /= s;
		}
		Mat3 Mat3::operator+(const Mat3& m) const
		{
			return Mat3(*this) += m;
		}
		Mat3 Mat3::operator-(const Mat3& m) const
		{
			return Mat3(*this) -= m;
		}
		Mat3 Mat3::operator*(const Mat3& m) const
		{
			Mat3 r;
			r.m_data[0] = m_data[0] * m.m_data[0] + m_data[3] * m.m_data[1] + m_data[6] * m.m_data[2];
			r.m_data[1] = m_data[1] * m.m_data[0] + m_data[4] * m.m_data[1] + m_data[7] * m.m_data[2];
			r.m_data[2] = m_data[2] * m.m_data[0] + m_data[5] * m.m_data[1] + m_data[8] * m.m_data[2];
			r.m_data[3] = m_data[0] * m.m_data[3] + m_data[3] * m.m_data[4] + m_data[6] * m.m_data[5];
			r.m_data[4] = m_data[1] * m.m_data[3] + m_data[4] * m.m_data[4] + m_data[7] * m.m_data[5];
			r.m_data[5] = m_data[2] * m.m_data[3] + m_data[5] * m.m_data[4] + m_data[8] * m.m_data[5];
			r.m_data[6] = m_data[0] * m.m_data[6] + m_data[3] * m.m_data[7] + m_data[6] * m.m_data[8];
			r.m_data[7] = m_data[1] * m.m_data[6] + m_data[4] * m.m_data[7] + m_data[7] * m.m_data[8];
			r.m_data[8] = m_data[2] * m.m_data[6] + m_data[5] * m.m_data[7] + m_data[8] * m.m_data[8];
			return r;
		}
		Vec3 Mat3::operator*(const Vec3& v) const
		{
			Vec3 r;
			r.m_data[0] = m_data[0] * v.m_data[0] + m_data[3] * v.m_data[1] + m_data[6] * v.m_data[2];
			r.m_data[1] = m_data[1] * v.m_data[0] + m_data[4] * v.m_data[1] + m_data[7] * v.m_data[2];
			r.m_data[2] = m_data[2] * v.m_data[0] + m_data[5] * v.m_data[1] + m_data[8] * v.m_data[2];
			return r;
		}

		Mat3& Mat3::operator+=(double s)
		{
			for(uint32_t i = 0; i < 9; i++)
			{
				m_data[i] += s;
			}
			return *this;
		}
		Mat3& Mat3::operator-=(double s)
		{
			for(uint32_t i = 0; i < 9; i++)
			{
				m_data[i] -= s;
			}
			return *this;
		}
		Mat3& Mat3::operator*=(double s)
		{
			for(uint32_t i = 0; i < 9; i++)
			{
				m_data[i] *= s;
			}
			return *this;
		}
		Mat3& Mat3::operator/=(double s)
		{
			for(uint32_t i = 0; i < 9; i++)
			{
				m_data[i] /= s;
			}
			return *this;
		}

		Mat3& Mat3::operator+=(const Mat3& m)
		{
			for(uint32_t i = 0; i < 9; i++)
			{
				m_data[i] += m.m_data[i];
			}
			return *this;
		}
		Mat3& Mat3::operator-=(const Mat3& m)
		{
			for(uint32_t i = 0; i < 9; i++)
			{
				m_data[i] -= m.m_data[i];
			}
			return *this;
		}

		double& Mat3::operator[](uint32_t i)
		{
			return m_data[i];
		}
		double& Mat3::operator()(uint32_t i, uint32_t j)
		{
			return m_data[i + 3 * j];
		}

		const double& Mat3::operator[](uint32_t i) const
		{
			return m_data[i];
		}
		const double& Mat3::operator()(uint32_t i, uint32_t j) const
		{
			return m_data[i + 3 * j];
		}

		//linear
		Mat3 Mat3::eye(void)
		{
			return Mat3({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
		}
		Vec3 Mat3::eigenvalues(void) const
		{
			Vec3 s;
			const double tl = lode_angle();
			const double I1 = invariant_I1();
			const double J2 = deviatoric_J2();
			s.m_data[2] = I1 / 3 + 2 * sqrt(J2 / 3) * cos(tl);
			s.m_data[0] = I1 / 3 + 2 * sqrt(J2 / 3) * cos(tl + 2 * M_PI / 3);
			s.m_data[1] = I1 / 3 + 2 * sqrt(J2 / 3) * cos(tl - 2 * M_PI / 3);
			return s;
		}
		Mat3 Mat3::inverse(void) const
		{
			Mat3 r;
			const double d = determinant();
			r.m_data[0] = (m_data[4] * m_data[8] - m_data[5] * m_data[7]) / d;
			r.m_data[1] = (m_data[2] * m_data[7] - m_data[1] * m_data[8]) / d;
			r.m_data[2] = (m_data[1] * m_data[5] - m_data[2] * m_data[4]) / d;
			r.m_data[3] = (m_data[5] * m_data[6] - m_data[3] * m_data[8]) / d;
			r.m_data[4] = (m_data[0] * m_data[8] - m_data[2] * m_data[6]) / d;
			r.m_data[5] = (m_data[2] * m_data[3] - m_data[0] * m_data[5]) / d;
			r.m_data[6] = (m_data[3] * m_data[7] - m_data[4] * m_data[6]) / d;
			r.m_data[7] = (m_data[6] * m_data[1] - m_data[0] * m_data[7]) / d;
			r.m_data[8] = (m_data[0] * m_data[4] - m_data[1] * m_data[3]) / d;
			return r;
		}
		Mat3 Mat3::transpose(void) const
		{
			Mat3 r(*this);
			r.m_data[1] = m_data[3];
			r.m_data[3] = m_data[1];
			r.m_data[2] = m_data[6];
			r.m_data[6] = m_data[2];
			r.m_data[5] = m_data[7];
			r.m_data[7] = m_data[5];
			return r;
		}

		double Mat3::trace(void) const
		{
			return m_data[0] + m_data[4] + m_data[8];
		}
		double Mat3::lode_angle(void) const
		{
			const double J2 = deviatoric_J2();
			const double J3 = deviatoric_J3();
			return fabs(J2) > 1e-5 ? acos(3 * sqrt(3) / 2 * J3 / sqrt(J2 * J2 * J2)) / 3 : 0;
		}
		double Mat3::determinant(void) const
		{
			return
				m_data[0] * (m_data[4] * m_data[8] - m_data[5] * m_data[7]) +
				m_data[3] * (m_data[2] * m_data[7] - m_data[1] * m_data[8]) +
				m_data[6] * (m_data[1] * m_data[5] - m_data[2] * m_data[4]);
		}
		double Mat3::invariant_I1(void) const
		{
			return m_data[0] + m_data[4] + m_data[8];
		}
		double Mat3::invariant_I2(void) const
		{
			return
				m_data[0] * m_data[4] + m_data[0] * m_data[8] + m_data[4] * m_data[8] -
				m_data[1] * m_data[1] - m_data[2] * m_data[2] - m_data[5] * m_data[5];
		}
		double Mat3::invariant_I3(void) const
		{
			return
				m_data[0] * m_data[4] * m_data[8] + 2 * m_data[1] * m_data[2] * m_data[5] -
				m_data[0] * m_data[5] * m_data[5] - m_data[4] * m_data[2] * m_data[2] - m_data[8] * m_data[1] * m_data[1];
		}
		double Mat3::deviatoric_J2(void) const
		{
			const double I1 = invariant_I1();
			const double I2 = invariant_I2();
			return I1 * I1 / 3 - I2;
		}
		double Mat3::deviatoric_J3(void) const
		{
			const double I1 = invariant_I1();
			const double I2 = invariant_I2();
			const double I3 = invariant_I3();
			return 2 * I1 * I1 * I1 / 27 - I1 * I2 / 3 + I3;
		}

		//friends
		Mat3 operator*(double s, const Mat3& m)
		{
			return Mat3(m) *= s;
		}
	}
}
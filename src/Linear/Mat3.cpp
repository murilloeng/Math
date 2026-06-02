//std
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Quat.hpp"
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Mat3.hpp"

namespace math
{
	//constructors
	Mat3::Mat3(mode mode) : Matrix(3, 3, mode)
	{
		return;
	}
	Mat3::Mat3(double* ptr) : Matrix(ptr, 3, 3)
	{
		return;
	}
	Mat3::Mat3(const Mat3& m) : Matrix(m)
	{
		return;
	}
	Mat3::Mat3(const double* ref) : Matrix(ref, 3, 3)
	{
		return;
	}
	Mat3::Mat3(std::initializer_list<double> list) : Matrix(3, 3)
	{
		if(list.size() != 9)
		{
			throw std::runtime_error("Mat3 constructor has incompatible dimensions!");
		}
		memcpy(m_data_ptr, std::data(list), list.size() * sizeof(double));
	}
	Mat3::Mat3(const Vec3& v1, const Vec3& v2, const Vec3& v3) : Matrix(3, 3)
	{
		memcpy(m_data_ptr + 0, v1.data(), 3 * sizeof(double));
		memcpy(m_data_ptr + 3, v2.data(), 3 * sizeof(double));
		memcpy(m_data_ptr + 6, v3.data(), 3 * sizeof(double));
	}
	Mat3::Mat3(const double* s1, const double* s2, const double* s3) : Matrix(3, 3)
	{
		memcpy(m_data_ptr + 0, s1, 3 * sizeof(double));
		memcpy(m_data_ptr + 3, s2, 3 * sizeof(double));
		memcpy(m_data_ptr + 6, s3, 3 * sizeof(double));
	}

	//destructor
	Mat3::~Mat3(void)
	{
		return;
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
		(Matrix&) r = (Matrix&) *this * m;
		return r;
	}
	Vec3 Mat3::operator*(const Vec3& v) const
	{
		Vec3 r;
		(Vector&) r = (Matrix&) * this * v;
		return r;
	}

	Mat3& Mat3::operator+=(double s)
	{
		(Matrix&) *this += s;
		return *this;
	}
	Mat3& Mat3::operator-=(double s)
	{
		(Matrix&) *this -= s;
		return *this;
	}
	Mat3& Mat3::operator*=(double s)
	{
		(Matrix&) *this *= s;
		return *this;
	}
	Mat3& Mat3::operator/=(double s)
	{
		(Matrix&) *this /= s;
		return *this;
	}

	Mat3& Mat3::operator=(const Mat3& m)
	{
		(Matrix&) *this = m;
		return *this;
	}

	Mat3& Mat3::operator+=(const Mat3& m)
	{
		(Matrix&) *this += m;
		return *this;
	}
	Mat3& Mat3::operator-=(const Mat3& m)
	{
		(Matrix&) *this -= m;
		return *this;
	}

	//linear
	Mat3 Mat3::eye(void)
	{
		Mat3 m;
		((Matrix&) m).eye();
		return m;
	}
	Vec3 Mat3::eigen(void) const
	{
		Vec3 s;
		const double tl = lode_angle();
		const double I1 = invariant_I1();
		const double J2 = deviatoric_J2();
		s[2] = I1 / 3 + 2 * sqrt(J2 / 3) * cos(tl);
		s[0] = I1 / 3 + 2 * sqrt(J2 / 3) * cos(tl + 2 * M_PI / 3);
		s[1] = I1 / 3 + 2 * sqrt(J2 / 3) * cos(tl - 2 * M_PI / 3);
		return s;
	}
	Mat3 Mat3::inverse(void) const
	{
		Mat3 r;
		const double d = determinant();
		r.m_data_ptr[0] = (m_data_ref[4] * m_data_ref[8] - m_data_ref[5] * m_data_ref[7]) / d;
		r.m_data_ptr[1] = (m_data_ref[2] * m_data_ref[7] - m_data_ref[1] * m_data_ref[8]) / d;
		r.m_data_ptr[2] = (m_data_ref[1] * m_data_ref[5] - m_data_ref[2] * m_data_ref[4]) / d;
		r.m_data_ptr[3] = (m_data_ref[5] * m_data_ref[6] - m_data_ref[3] * m_data_ref[8]) / d;
		r.m_data_ptr[4] = (m_data_ref[0] * m_data_ref[8] - m_data_ref[2] * m_data_ref[6]) / d;
		r.m_data_ptr[5] = (m_data_ref[2] * m_data_ref[3] - m_data_ref[0] * m_data_ref[5]) / d;
		r.m_data_ptr[6] = (m_data_ref[3] * m_data_ref[7] - m_data_ref[4] * m_data_ref[6]) / d;
		r.m_data_ptr[7] = (m_data_ref[6] * m_data_ref[1] - m_data_ref[0] * m_data_ref[7]) / d;
		r.m_data_ptr[8] = (m_data_ref[0] * m_data_ref[4] - m_data_ref[1] * m_data_ref[3]) / d;
		return r;
	}
	Mat3 Mat3::transpose(void) const
	{
		Mat3 r;
		(Matrix&) r = ((Matrix&) *this).transpose();
		return r;
	}

	double Mat3::trace(void) const
	{
		return m_data_ref[0] + m_data_ref[4] + m_data_ref[8];
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
			m_data_ref[0] * (m_data_ref[4] * m_data_ref[8] - m_data_ref[5] * m_data_ref[7]) +
			m_data_ref[3] * (m_data_ref[2] * m_data_ref[7] - m_data_ref[1] * m_data_ref[8]) +
			m_data_ref[6] * (m_data_ref[1] * m_data_ref[5] - m_data_ref[2] * m_data_ref[4]);
	}
	double Mat3::invariant_I1(void) const
	{
		return m_data_ref[0] + m_data_ref[4] + m_data_ref[8];
	}
	double Mat3::invariant_I2(void) const
	{
		return
			m_data_ref[0] * m_data_ref[4] + m_data_ref[0] * m_data_ref[8] + m_data_ref[4] * m_data_ref[8] -
			m_data_ref[1] * m_data_ref[1] - m_data_ref[2] * m_data_ref[2] - m_data_ref[5] * m_data_ref[5];
	}
	double Mat3::invariant_I3(void) const
	{
		return
			m_data_ref[0] * m_data_ref[4] * m_data_ref[8] - 2 * m_data_ref[1] * m_data_ref[2] * m_data_ref[5] -
			m_data_ref[0] * m_data_ref[5] * m_data_ref[5] - m_data_ref[4] * m_data_ref[2] * m_data_ref[2] - m_data_ref[8] * m_data_ref[1] * m_data_ref[1];
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

	//rotation
	Vec3 Mat3::rotation(void) const
	{
		return quaternion().pseudo();
	}
	Quat Mat3::quaternion(void) const
	{
		Quat q;
		const double d0 = m_data_ref[0];
		const double d1 = m_data_ref[4];
		const double d2 = m_data_ref[8];
		const double tr = d0 + d1 + d2;
		if(tr >= d0 && tr >= d1 && tr >= d2)
		{
			q[0] = sqrt(1 + tr) / 2;
			q[1] = (m_data_ref[2 + 3 * 1] - m_data_ref[1 + 3 * 2]) / (4 * q[0]);
			q[2] = (m_data_ref[0 + 3 * 2] - m_data_ref[2 + 3 * 0]) / (4 * q[0]);
			q[3] = (m_data_ref[1 + 3 * 0] - m_data_ref[0 + 3 * 1]) / (4 * q[0]);
		}
		else if(d0 >= d1 && d0 >= d2)
		{
			q[1] = sqrt(1 - tr + 2 * d0) / 2;
			q[0] = (m_data_ref[2 + 3 * 1] - m_data_ref[1 + 3 * 2]) / (4 * q[1]);
			q[2] = (m_data_ref[0 + 3 * 1] + m_data_ref[1 + 3 * 0]) / (4 * q[1]);
			q[3] = (m_data_ref[0 + 3 * 2] + m_data_ref[2 + 3 * 0]) / (4 * q[1]);
		}
		else if(d1 >= d2)
		{
			q[2] = sqrt(1 - tr + 2 * d1) / 2;
			q[0] = (m_data_ref[0 + 3 * 2] - m_data_ref[2 + 3 * 0]) / (4 * q[2]);
			q[1] = (m_data_ref[0 + 3 * 1] + m_data_ref[1 + 3 * 0]) / (4 * q[2]);
			q[3] = (m_data_ref[1 + 3 * 2] + m_data_ref[2 + 3 * 1]) / (4 * q[2]);
		}
		else
		{
			q[3] = sqrt(1 - tr + 2 * d2) / 2;
			q[0] = (m_data_ref[1 + 3 * 0] - m_data_ref[0 + 3 * 1]) / (4 * q[3]);
			q[1] = (m_data_ref[0 + 3 * 2] + m_data_ref[2 + 3 * 0]) / (4 * q[3]);
			q[2] = (m_data_ref[1 + 3 * 2] + m_data_ref[2 + 3 * 1]) / (4 * q[3]);
		}
		return q;
	}

	//friends
	Mat3 operator*(double s, const Mat3& m)
	{
		return Mat3(m) *= s;
	}
}
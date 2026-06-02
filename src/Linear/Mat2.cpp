//std
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Vec2.hpp"
#include "Math/inc/Linear/Mat2.hpp"

namespace math
{
	//constructors
	Mat2::Mat2(mode init) : Matrix(2, 2, init)
	{
		return;
	}
	Mat2::Mat2(double* ptr) : Matrix(ptr, 2, 2)
	{
		return;
	}
	Mat2::Mat2(const Mat2& m) : Matrix(m)
	{
		return;
	}
	Mat2::Mat2(const double* ref) : Matrix(ref, 2, 2)
	{
		return;
	}
	Mat2::Mat2(const Vec2& v1, const Vec2& v2) : Matrix(2, 2)
	{
		memcpy(m_data_ptr + 0, v1.data(), 2 * sizeof(double));
		memcpy(m_data_ptr + 2, v2.data(), 2 * sizeof(double));
	}
	Mat2::Mat2(const double* s1, const double* s2) : Matrix(2, 2)
	{
		memcpy(m_data_ptr + 0, s1, 2 * sizeof(double));
		memcpy(m_data_ptr + 2, s2, 2 * sizeof(double));
	}
	Mat2::Mat2(std::initializer_list<double> list) : Matrix(2, 2)
	{
		if(list.size() != 4)
		{
			throw std::runtime_error("Mat2 constructor has incompatible dimensions!");
		}
		memcpy(m_data_ptr, std::data(list), list.size() * sizeof(double));
	}

	//destructor
	Mat2::~Mat2(void)
	{
		return;
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
		r.m_data_ptr[0] = m_data_ref[0] * m.m_data_ref[0] + m_data_ref[2] * m.m_data_ref[1];
		r.m_data_ptr[1] = m_data_ref[1] * m.m_data_ref[0] + m_data_ref[3] * m.m_data_ref[1];
		r.m_data_ptr[2] = m_data_ref[0] * m.m_data_ref[2] + m_data_ref[2] * m.m_data_ref[3];
		r.m_data_ptr[3] = m_data_ref[1] * m.m_data_ref[2] + m_data_ref[3] * m.m_data_ref[3];
		return r;
	}
	Vec2 Mat2::operator*(const Vec2& v) const
	{
		Vec2 r;
		r[0] = m_data_ref[0] * v[0] + m_data_ref[2] * v[1];
		r[1] = m_data_ref[1] * v[0] + m_data_ref[3] * v[1];
		return r;
	}

	Mat2& Mat2::operator+=(double s)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_data_ptr[i] += s;
		}
		return *this;
	}
	Mat2& Mat2::operator-=(double s)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_data_ptr[i] -= s;
		}
		return *this;
	}
	Mat2& Mat2::operator*=(double s)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_data_ptr[i] *= s;
		}
		return *this;
	}
	Mat2& Mat2::operator/=(double s)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_data_ptr[i] /= s;
		}
		return *this;
	}

	Mat2& Mat2::operator=(const Mat2& m)
	{
		memcpy(m_data_ptr, m.m_data_ref, 4 * sizeof(double));
		return *this;
	}

	Mat2& Mat2::operator+=(const Mat2& m)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_data_ptr[i] += m.m_data_ref[i];
		}
		return *this;
	}
	Mat2& Mat2::operator-=(const Mat2& m)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_data_ptr[i] -= m.m_data_ref[i];
		}
		return *this;
	}

	double& Mat2::operator[](uint32_t i)
	{
		return m_data_ptr[i];
	}
	double& Mat2::operator()(uint32_t i)
	{
		return m_data_ptr[i];
	}
	double& Mat2::operator()(uint32_t i, uint32_t j)
	{
		return m_data_ptr[i + 2 * j];
	}

	const double& Mat2::operator[](uint32_t i) const
	{
		return m_data_ref[i];
	}
	const double& Mat2::operator()(uint32_t i) const
	{
		return m_data_ref[i];
	}
	const double& Mat2::operator()(uint32_t i, uint32_t j) const
	{
		return m_data_ref[i + 2 * j];
	}

	//linear
	Mat2 Mat2::eye(void)
	{
		Mat2 m;
		m.Matrix::eye();
		return m;
	}
	Vec2 Mat2::eigen(void) const
	{
		Vec2 s;
		const double a = m_data_ref[0] + m_data_ref[3];
		const double b = m_data_ref[0] - m_data_ref[3];
		const double c = m_data_ref[1] * m_data_ref[2];
		s[0] = a / 2 - sqrt(b * b + 4 * c) / 2;
		s[1] = a / 2 + sqrt(b * b + 4 * c) / 2;
		return s;
	}
	Mat2 Mat2::inverse(void) const
	{
		Mat2 r;
		r.m_data_ptr[0] = +m_data_ref[3] / (m_data_ref[0] * m_data_ref[3] - m_data_ref[1] * m_data_ref[2]);
		r.m_data_ptr[1] = -m_data_ref[1] / (m_data_ref[0] * m_data_ref[3] - m_data_ref[1] * m_data_ref[2]);
		r.m_data_ptr[2] = -m_data_ref[2] / (m_data_ref[0] * m_data_ref[3] - m_data_ref[1] * m_data_ref[2]);
		r.m_data_ptr[3] = +m_data_ref[0] / (m_data_ref[0] * m_data_ref[3] - m_data_ref[1] * m_data_ref[2]);
		return r;
	}
	Mat2 Mat2::transpose(void) const
	{
		Mat2 r;
		r.m_data_ptr[0] = m_data_ref[0];
		r.m_data_ptr[1] = m_data_ref[2];
		r.m_data_ptr[2] = m_data_ref[1];
		r.m_data_ptr[3] = m_data_ref[3];
		return r;
	}

	double Mat2::trace(void) const
	{
		return m_data_ref[0] + m_data_ref[3];
	}
	double Mat2::determinant(void) const
	{
		return m_data_ref[0] * m_data_ref[3] - m_data_ref[1] * m_data_ref[2];
	}
	double Mat2::invariant_I1(void) const
	{
		return m_data_ref[0] + m_data_ref[3];
	}
	double Mat2::invariant_I2(void) const
	{
		return m_data_ref[0] * m_data_ref[3] - m_data_ref[1] * m_data_ref[2];
	}
	double Mat2::deviatoric_J2(void) const
	{
		return -(m_data_ref[0] - m_data_ref[3]) * (m_data_ref[0] - m_data_ref[3]) / 4 - m_data_ref[1] * m_data_ref[2];
	}

	//friends
	Mat2 operator*(double s, const Mat2& m)
	{
		return Mat2(m) *= s;
	}
}
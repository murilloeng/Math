//std
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Mat3.hpp"
#include "Math/inc/Linear/Mat4.hpp"

namespace math
{
	//constructors
	Mat4::Mat4(mode mode) : Matrix(4, 4, mode)
	{
		return;
	}
	Mat4::Mat4(double* ptr) : Matrix(ptr, 4, 4)
	{
		return;
	}
	Mat4::Mat4(const Mat4& m) : Matrix(m)
	{
		return;
	}
	Mat4::Mat4(const double* ref) : Matrix(ref, 4, 4)
	{
		return;
	}
	Mat4::Mat4(const Mat3& M, const Vec3& v) : Matrix(4, 4)
	{
		m_data_ptr[3 + 4 * 0] = 0;
		m_data_ptr[3 + 4 * 1] = 0;
		m_data_ptr[3 + 4 * 2] = 0;
		m_data_ptr[3 + 4 * 3] = 1;
		for(uint32_t i = 0; i < 3; i++)
		{
			m_data_ptr[i + 4 * 3] = v[i];
			for(uint32_t j = 0; j < 3; j++)
			{
				m_data_ptr[i + 4 * j] = M[i + 3 * j];
			}
		}
	}
	Mat4::Mat4(std::initializer_list<double> list) : Matrix(4, 4)
	{
		if(list.size() != 16)
		{
			throw std::runtime_error("Mat4 constructor has incompatible dimensions!");
		}
		memcpy(m_data_ptr, std::data(list), list.size() * sizeof(double));
	}

	//destructor
	Mat4::~Mat4(void)
	{
		return;
	}

	//operators
	Mat4 Mat4::operator+(void) const
	{
		return *this;
	}
	Mat4 Mat4::operator-(void) const
	{
		return Mat4(*this) *= -1;
	}

	Mat4 Mat4::operator/(double s) const
	{
		return Mat4(*this) /= s;
	}
	Mat4 Mat4::operator+(const Mat4& m) const
	{
		return Mat4(*this) += m;
	}
	Mat4 Mat4::operator-(const Mat4& m) const
	{
		return Mat4(*this) -= m;
	}
	Mat4 Mat4::operator*(const Mat4& m) const
	{
		Mat4 r;
		((Matrix&) r) = ((Matrix&) *this) * m;
		return r;
	}
	Vec3 Mat4::operator*(const Vec3& v) const
	{
		Vec3 r;
		for(uint32_t i = 0; i < 3; i++)
		{
			r[i] = m_data_ref[i + 4 * 3];
			for(uint32_t j = 0; j < 3; j++)
			{
				r[i] += m_data_ref[i + 4 * j] * v[j];
			}
		}
		return r;
	}

	Mat4& Mat4::operator+=(double s)
	{
		(Matrix&) *this += s;
		return *this;
	}
	Mat4& Mat4::operator-=(double s)
	{
		(Matrix&) *this -= s;
		return *this;
	}
	Mat4& Mat4::operator*=(double s)
	{
		(Matrix&) *this *= s;
		return *this;
	}
	Mat4& Mat4::operator/=(double s)
	{
		(Matrix&) *this /= s;
		return *this;
	}

	Mat4& Mat4::operator=(const Mat4& m)
	{
		(Matrix&) *this = m;
		return *this;
	}

	Mat4& Mat4::operator+=(const Mat4& m)
	{
		(Matrix&) *this += m;
		return *this;
	}
	Mat4& Mat4::operator-=(const Mat4& m)
	{
		(Matrix&) *this -= m;
		return *this;
	}

	//linear
	Mat4 Mat4::eye(void)
	{
		Mat4 m;
		((Matrix&) m).eye();
		return m;
	}
	Mat4 Mat4::transpose(void) const
	{
		Mat4 r;
		(Matrix&) r = ((Matrix&) *this).transpose();
		return r;
	}

	//friends
	Mat4 operator*(double s, const Mat4& m)
	{
		return Mat4(m) *= s;
	}
}
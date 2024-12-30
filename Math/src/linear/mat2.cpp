//std
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//math
#include "Math/Math/inc/linear/vec2.hpp"
#include "Math/Math/inc/linear/mat2.hpp"

namespace math
{
	//constructors
	mat2::mat2(mode init) : matrix(2, 2, init)
	{
		return;
	}
	mat2::mat2(double* ptr) : matrix(ptr, 2, 2)
	{
		return;
	}
	mat2::mat2(const mat2& m) : matrix(m)
	{
		return;
	}
	mat2::mat2(const double* ref) : matrix(ref, 2, 2)
	{
		return;
	}
	mat2::mat2(const vec2& v1, const vec2& v2) : matrix(2, 2)
	{
		memcpy(m_ptr + 0, v1.data(), 2 * sizeof(double));
		memcpy(m_ptr + 2, v2.data(), 2 * sizeof(double));
	}
	mat2::mat2(const double* s1, const double* s2) : matrix(2, 2)
	{
		memcpy(m_ptr + 0, s1, 2 * sizeof(double));
		memcpy(m_ptr + 2, s2, 2 * sizeof(double));
	}
	mat2::mat2(std::initializer_list<double> list) : matrix(2, 2)
	{
		if(list.size() != 4)
		{
			fprintf(stderr, "Error: mat2 constructor with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		memcpy(m_ptr, std::data(list), list.size() * sizeof(double));
	}

	//destructor
	mat2::~mat2(void)
	{
		return;
	}

	//operators
	mat2 mat2::operator+(void) const
	{
		return *this;
	}
	mat2 mat2::operator-(void) const
	{
		return mat2(*this) *= -1;
	}

	mat2 mat2::operator/(double s) const
	{
		return mat2(*this) /= s;
	}
	mat2 mat2::operator+(const mat2& m) const
	{
		return mat2(*this) += m;
	}
	mat2 mat2::operator-(const mat2& m) const
	{
		return mat2(*this) -= m;
	}
	mat2 mat2::operator*(const mat2& m) const
	{
		mat2 r;
		r.m_ptr[0] = m_ref[0] * m.m_ref[0] + m_ref[2] * m.m_ref[1];
		r.m_ptr[1] = m_ref[1] * m.m_ref[0] + m_ref[3] * m.m_ref[1];
		r.m_ptr[2] = m_ref[0] * m.m_ref[2] + m_ref[2] * m.m_ref[3];
		r.m_ptr[3] = m_ref[1] * m.m_ref[2] + m_ref[3] * m.m_ref[3];
		return r;
	}
	vec2 mat2::operator*(const vec2& v) const
	{
		vec2 r;
		r[0] = m_ref[0] * v[0] + m_ref[2] * v[1];
		r[1] = m_ref[1] * v[0] + m_ref[3] * v[1];
		return r;
	}

	mat2& mat2::operator+=(double s)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_ptr[i] += s;
		}
		return *this;
	}
	mat2& mat2::operator-=(double s)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_ptr[i] -= s;
		}
		return *this;
	}
	mat2& mat2::operator*=(double s)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_ptr[i] *= s;
		}
		return *this;
	}
	mat2& mat2::operator/=(double s)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_ptr[i] /= s;
		}
		return *this;
	}

	mat2& mat2::operator=(const mat2& m)
	{
		memcpy(m_ptr, m.m_ref, 4 * sizeof(double));
		return *this;
	}

	mat2& mat2::operator+=(const mat2& m)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_ptr[i] += m.m_ref[i];
		}
		return *this;
	}
	mat2& mat2::operator-=(const mat2& m)
	{
		for(uint32_t i = 0; i < 4; i++)
		{
			m_ptr[i] -= m.m_ref[i];
		}
		return *this;
	}

	double& mat2::operator[](uint32_t i)
	{
		return m_ptr[i];
	}
	double& mat2::operator()(uint32_t i)
	{
		return m_ptr[i];
	}
	double& mat2::operator()(uint32_t i, uint32_t j)
	{
		return m_ptr[i + 2 * j];
	}

	const double& mat2::operator[](uint32_t i) const
	{
		return m_ref[i];
	}
	const double& mat2::operator()(uint32_t i) const
	{
		return m_ref[i];
	}
	const double& mat2::operator()(uint32_t i, uint32_t j) const
	{
		return m_ref[i + 2 * j];
	}

	//linear
	mat2 mat2::eye(void)
	{
		mat2 m;
		m.matrix::eye();
		return m;
	}
	vec2 mat2::eigen(void) const
	{
		vec2 s;
		const double a = m_ref[0] + m_ref[3];
		const double b = m_ref[0] - m_ref[3];
		const double c = m_ref[1] * m_ref[2];
		s[0] = a / 2 - sqrt(b * b + 4 * c) / 2;
		s[1] = a / 2 + sqrt(b * b + 4 * c) / 2;
		return s;
	}
	mat2 mat2::inverse(void) const
	{
		mat2 r;
		r.m_ptr[0] = +m_ref[3] / (m_ref[0] * m_ref[3] - m_ref[1] * m_ref[2]);
		r.m_ptr[1] = -m_ref[1] / (m_ref[0] * m_ref[3] - m_ref[1] * m_ref[2]);
		r.m_ptr[2] = -m_ref[2] / (m_ref[0] * m_ref[3] - m_ref[1] * m_ref[2]);
		r.m_ptr[3] = +m_ref[0] / (m_ref[0] * m_ref[3] - m_ref[1] * m_ref[2]);
		return r;
	}
	mat2 mat2::transpose(void) const
	{
		mat2 r;
		r.m_ptr[0] = m_ref[0];
		r.m_ptr[1] = m_ref[2];
		r.m_ptr[2] = m_ref[1];
		r.m_ptr[3] = m_ref[3];
		return r;
	}

	double mat2::trace(void) const
	{
		return m_ref[0] + m_ref[3];
	}
	double mat2::determinant(void) const
	{
		return m_ref[0] * m_ref[3] - m_ref[1] * m_ref[2];
	}
	double mat2::invariant_I1(void) const
	{
		return m_ref[0] + m_ref[3];
	}
	double mat2::invariant_I2(void) const
	{
		return m_ref[0] * m_ref[3] - m_ref[1] * m_ref[2];
	}
	double mat2::deviatoric_J2(void) const
	{
		return -(m_ref[0] - m_ref[3]) * (m_ref[0] - m_ref[3]) / 4 - m_ref[1] * m_ref[2];
	}

	//friends
	mat2 operator*(double s, const mat2& m)
	{
		return mat2(m) *= s;
	}
}
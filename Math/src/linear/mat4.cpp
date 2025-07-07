//std
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//math
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat4.hpp"

namespace math
{
	//constructors
	mat4::mat4(mode mode) : matrix(4, 4, mode)
	{
		return;
	}
	mat4::mat4(double* ptr) : matrix(ptr, 4, 4)
	{
		return;
	}
	mat4::mat4(const mat4& m) : matrix(m)
	{
		return;
	}
	mat4::mat4(const double* ref) : matrix(ref, 4, 4)
	{
		return;
	}
	mat4::mat4(std::initializer_list<double> list) : matrix(4, 4)
	{
		if(list.size() != 16)
		{
			fprintf(stderr, "Error: mat4 constructor with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		memcpy(m_data_ptr, std::data(list), list.size() * sizeof(double));
	}

	//destructor
	mat4::~mat4(void)
	{
		return;
	}

	//operators
	mat4 mat4::operator+(void) const
	{
		return *this;
	}
	mat4 mat4::operator-(void) const
	{
		return mat4(*this) *= -1;
	}

	mat4 mat4::operator/(double s) const
	{
		return mat4(*this) /= s;
	}
	mat4 mat4::operator+(const mat4& m) const
	{
		return mat4(*this) += m;
	}
	mat4 mat4::operator-(const mat4& m) const
	{
		return mat4(*this) -= m;
	}
	mat4 mat4::operator*(const mat4& m) const
	{
		mat4 r;
		((matrix&) r) = ((matrix&) *this) * m;
		return r;
	}
	vec3 mat4::operator*(const vec3& v) const
	{
		vec3 r;
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

	mat4& mat4::operator+=(double s)
	{
		(matrix&) *this += s;
		return *this;
	}
	mat4& mat4::operator-=(double s)
	{
		(matrix&) *this -= s;
		return *this;
	}
	mat4& mat4::operator*=(double s)
	{
		(matrix&) *this *= s;
		return *this;
	}
	mat4& mat4::operator/=(double s)
	{
		(matrix&) *this /= s;
		return *this;
	}

	mat4& mat4::operator=(const mat4& m)
	{
		(matrix&) *this = m;
		return *this;
	}

	mat4& mat4::operator+=(const mat4& m)
	{
		(matrix&) *this += m;
		return *this;
	}
	mat4& mat4::operator-=(const mat4& m)
	{
		(matrix&) *this -= m;
		return *this;
	}

	//linear
	mat4 mat4::eye(void)
	{
		mat4 m;
		((matrix&) m).eye();
		return m;
	}
	mat4 mat4::transpose(void) const
	{
		mat4 r;
		(matrix&) r = ((matrix&) *this).transpose();
		return r;
	}

	//friends
	mat4 operator*(double s, const mat4& m)
	{
		return mat4(m) *= s;
	}
}
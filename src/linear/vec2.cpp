//std
#include <cmath>

//math
#include "Math/inc/linear/vec2.hpp"

namespace math
{
	//constructors
	vec2::vec2(void) : vector(2)
	{
		return;
	}
	vec2::vec2(double* ptr) : vector(ptr, 2)
	{
		return;
	}
	vec2::vec2(const vec2& v) : vector(v)
	{
		return;
	}
	vec2::vec2(const double* ref) : vector(ref, 2)
	{
		return;
	}
	vec2::vec2(double v0, double v1) : vector({v0, v1})
	{
		return;
	}

	//destructor
	vec2::~vec2(void)
	{
		return;
	}

	//operators
	vec2 vec2::operator+(void) const
	{
		return *this;
	}
	vec2 vec2::operator-(void) const
	{
		return vec2(*this) *= -1;
	}

	vec2 vec2::operator/(double s) const
	{
		return vec2(*this) /= s;
	}
	vec2 vec2::operator+(const vec2& v) const
	{
		return vec2(*this) += v;
	}
	vec2 vec2::operator-(const vec2& v) const
	{
		return vec2(*this) -= v;
	}

	vec2& vec2::operator+=(double s)
	{
		m_ptr[0] += s;
		m_ptr[1] += s;
		return *this;
	}
	vec2& vec2::operator-=(double s)
	{
		m_ptr[0] -= s;
		m_ptr[1] -= s;
		return *this;
	}
	vec2& vec2::operator*=(double s)
	{
		m_ptr[0] *= s;
		m_ptr[1] *= s;
		return *this;
	}
	vec2& vec2::operator/=(double s)
	{
		m_ptr[0] /= s;
		m_ptr[1] /= s;
		return *this;
	}

	vec2& vec2::operator=(const vec2& v)
	{
		m_ptr[0] = v.m_ref[0];
		m_ptr[1] = v.m_ref[1];
		return *this;
	}

	vec2& vec2::operator+=(const vec2& v)
	{
		m_ptr[0] += v.m_ref[0];
		m_ptr[1] += v.m_ref[1];
		return *this;
	}
	vec2& vec2::operator-=(const vec2& v)
	{
		m_ptr[0] -= v.m_ref[0];
		m_ptr[1] -= v.m_ref[1];
		return *this;
	}

	double& vec2::operator[](uint32_t i)
	{
		return m_ptr[i];
	}
	double& vec2::operator()(uint32_t i)
	{
		return m_ptr[i];
	}

	const double& vec2::operator[](uint32_t i) const
	{
		return m_ref[i];
	}
	const double& vec2::operator()(uint32_t i) const
	{
		return m_ref[i];
	}

	//linear
	vec2& vec2::normalize(void)
	{
		return *this /= norm();
	}
	vec2& vec2::project(const vec2& v)
	{
		return *this -= inner(v) * v;
	}

	double vec2::inner(const vec2& v) const
	{
		return m_ref[0] * v.m_ref[0] + m_ref[1] * v.m_ref[1];
	}
	double vec2::cross(const vec2& v) const
	{
		return m_ref[0] * v.m_ref[1] - m_ref[1] * v.m_ref[0];
	}

	//friends
	vec2 operator*(double s, const vec2& v)
	{
		return vec2(v) *= s;
	}
}
//std
#include <cmath>

//Math
#include "Math/inc/Linear/Vec2.hpp"

namespace math
{
	//constructors
	Vec2::Vec2(void) : Vector(2)
	{
		return;
	}
	Vec2::Vec2(double* ptr) : Vector(ptr, 2)
	{
		return;
	}
	Vec2::Vec2(const Vec2& v) : Vector(v)
	{
		return;
	}
	Vec2::Vec2(const double* ref) : Vector(ref, 2)
	{
		return;
	}
	Vec2::Vec2(double v0, double v1) : Vector({v0, v1})
	{
		return;
	}

	//destructor
	Vec2::~Vec2(void)
	{
		return;
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
		m_data_ptr[0] += s;
		m_data_ptr[1] += s;
		return *this;
	}
	Vec2& Vec2::operator-=(double s)
	{
		m_data_ptr[0] -= s;
		m_data_ptr[1] -= s;
		return *this;
	}
	Vec2& Vec2::operator*=(double s)
	{
		m_data_ptr[0] *= s;
		m_data_ptr[1] *= s;
		return *this;
	}
	Vec2& Vec2::operator/=(double s)
	{
		m_data_ptr[0] /= s;
		m_data_ptr[1] /= s;
		return *this;
	}

	Vec2& Vec2::operator=(const Vec2& v)
	{
		m_data_ptr[0] = v.m_data_ref[0];
		m_data_ptr[1] = v.m_data_ref[1];
		return *this;
	}

	Vec2& Vec2::operator+=(const Vec2& v)
	{
		m_data_ptr[0] += v.m_data_ref[0];
		m_data_ptr[1] += v.m_data_ref[1];
		return *this;
	}
	Vec2& Vec2::operator-=(const Vec2& v)
	{
		m_data_ptr[0] -= v.m_data_ref[0];
		m_data_ptr[1] -= v.m_data_ref[1];
		return *this;
	}

	double& Vec2::operator[](uint32_t i)
	{
		return m_data_ptr[i];
	}
	double& Vec2::operator()(uint32_t i)
	{
		return m_data_ptr[i];
	}

	const double& Vec2::operator[](uint32_t i) const
	{
		return m_data_ref[i];
	}
	const double& Vec2::operator()(uint32_t i) const
	{
		return m_data_ref[i];
	}

	//linear
	Vec2& Vec2::normalize(void)
	{
		return *this /= norm();
	}
	Vec2& Vec2::project(const Vec2& v)
	{
		return *this -= inner(v) * v;
	}

	double Vec2::inner(const Vec2& v) const
	{
		return m_data_ref[0] * v.m_data_ref[0] + m_data_ref[1] * v.m_data_ref[1];
	}
	double Vec2::cross(const Vec2& v) const
	{
		return m_data_ref[0] * v.m_data_ref[1] - m_data_ref[1] * v.m_data_ref[0];
	}

	//friends
	Vec2 operator*(double s, const Vec2& v)
	{
		return Vec2(v) *= s;
	}
}
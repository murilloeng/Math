//std
#include <cmath>

//Math
#include "Math/inc/Groups/ASO2.hpp"
#include "Math/inc/Groups/GSO2.hpp"

namespace math
{
	//constructor
	GSO2::GSO2(double data) : m_data{data}
	{
		return;
	}

	//actions
	Matrix GSO2::matrix(void) const
	{
		const double c = cos(m_data);
		const double s = sin(m_data);
		return Matrix({{c, -s}, {s, c}});
	}
	Vector GSO2::action(const Vector& v) const
	{
		const double c = cos(m_data);
		const double s = sin(m_data);
		return Vector({c * v[0] - s * v[1], s * v[0] + c * v[1]});
	}

	//operations
	ASO2 GSO2::log(void) const
	{
		return ASO2(m_data);
	}
	GSO2 GSO2::inverse(void) const
	{
		return GSO2(-m_data);
	}
	GSO2 GSO2::operator*(const GSO2& element) const
	{
		return GSO2(m_data + element.m_data);
	}
	GSO2 GSO2::interpolate(const GSO2& element, double t) const
	{
		return GSO2((1 - t) * m_data + t * element.m_data);
	}
}
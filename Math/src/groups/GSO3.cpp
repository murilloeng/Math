//std
#include <cmath>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/groups/ASO3.hpp"
#include "Math/Math/inc/groups/GSO3.hpp"

namespace math
{
	namespace groups
	{
		//constructor
		GSO3::GSO3(void)
		{
			return;
		}
		GSO3::GSO3(quat quaternion) : m_quaternion{quaternion}
		{
			return;
		}
		
		//destructor
		GSO3::~GSO3(void)
		{
			return;
		}

		//inverse
		GSO3 GSO3::inverse(void) const
		{
			return m_quaternion.conjugate();
		}

		//logarithm
		ASO3 GSO3::logarithm(void) const
		{
			//angle
			const double t = 2 * acos(bound(m_quaternion[0]));
			//vector
			vec3 algebra;
			const double s = sin(t / 2);
			algebra[0] = s ? t * m_quaternion[1] / s : 0;
			algebra[1] = s ? t * m_quaternion[2] / s : 0;
			algebra[2] = s ? t * m_quaternion[3] / s : 0;
			//return
			return algebra;
		}

		//matrix
		mat3 GSO3::matrix_form(void) const
		{
			mat3 matrix;
			const double* p = m_quaternion.data();
			matrix[1 + 3 * 0] = 2 * (p[1] * p[2] + p[0] * p[3]);
			matrix[2 + 3 * 0] = 2 * (p[1] * p[3] - p[0] * p[2]);
			matrix[0 + 3 * 1] = 2 * (p[1] * p[2] - p[0] * p[3]);
			matrix[2 + 3 * 1] = 2 * (p[2] * p[3] + p[0] * p[1]);
			matrix[0 + 3 * 2] = 2 * (p[1] * p[3] + p[0] * p[2]);
			matrix[1 + 3 * 2] = 2 * (p[2] * p[3] - p[0] * p[1]);
			matrix[0 + 3 * 0] = p[0] * p[0] + p[1] * p[1] - p[2] * p[2] - p[3] * p[3];
			matrix[1 + 3 * 1] = p[0] * p[0] - p[1] * p[1] + p[2] * p[2] - p[3] * p[3];
			matrix[2 + 3 * 2] = p[0] * p[0] - p[1] * p[1] - p[2] * p[2] + p[3] * p[3];
			return matrix;
		}

		//quaternion
		quat& GSO3::quaternion(void)
		{
			return m_quaternion;
		}
		const quat& GSO3::quaternion(void) const
		{
			return m_quaternion;
		}

		//operators
		vec3 GSO3::operator*(const vec3& v) const
		{
			//data
			vec3 r;
			const double s = m_quaternion[0];
			const vec3 x(m_quaternion.data() + 1);
			//vector
			const double b = 2 * x.inner(v);
			const double a = s * s - x.inner(x);
			r[0] = a * v[0] + b * x[0] + 2 * s * (x[1] * v[2] - x[2] * v[1]);
			r[1] = a * v[1] + b * x[1] + 2 * s * (x[2] * v[0] - x[0] * v[2]);
			r[2] = a * v[2] + b * x[2] + 2 * s * (x[0] * v[1] - x[1] * v[0]);
			//return
			return r;
		}
		GSO3 GSO3::operator*(const GSO3& group) const
		{
			//data
			GSO3 r;
			//product
			r.m_quaternion = m_quaternion * group.m_quaternion;
			//return
			return r;
		}
	}
}
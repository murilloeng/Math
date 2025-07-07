//math
#include "Math/Math/inc/linear/mat3.hpp"
#include "Math/Math/inc/groups/ASO3.hpp"
#include "Math/Math/inc/groups/ASE3.hpp"
#include "Math/Math/inc/groups/GSO3.hpp"
#include "Math/Math/inc/groups/GSE3.hpp"

namespace math
{
	namespace groups
	{
		//constructors
		GSE3::GSE3(void)
		{
			return;
		}
		GSE3::GSE3(vec3 vector, quat quaternion) : m_vector{vector}, m_quaternion{quaternion}
		{
			return;
		}

		//destructor
		GSE3::~GSE3(void)
		{
			return;
		}

		//inverse
		GSE3 GSE3::inverse(void) const
		{
			return GSE3(-m_quaternion.conjugate(m_vector), m_quaternion.conjugate());
		}

		//logarithm
		ASE3 GSE3::logarithm(void) const
		{
			ASE3 object;
			object.m_vector_w = GSO3(m_quaternion).logarithm().vector();
			object.m_vector_u = ASO3(object.m_vector_w).tangent_inverse().transpose() * m_vector;
			return object;
		}

		//matrix
		mat4 GSE3::matrix_form(void) const
		{
			return mat4(GSO3(m_quaternion).matrix_form(), m_vector);
		}

		//vector
		vec3& GSE3::vector(void)
		{
			return m_vector;
		}
		const vec3& GSE3::vector(void) const
		{
			return m_vector;
		}

		//quaternion
		quat& GSE3::quaternion(void)
		{
			return m_quaternion;
		}
		const quat& GSE3::quaternion(void) const
		{
			return m_quaternion;
		}

		//operators
		vec3 GSE3::operator*(const vec3& v) const
		{
			return GSO3(m_quaternion) * v + m_vector;
		}
		GSE3 GSE3::operator*(const GSE3& object) const
		{
			const quat quaternion = m_quaternion * object.m_quaternion;
			const vec3 vector = m_vector + GSO3(m_quaternion) * object.m_vector;
			return GSE3(vector, quaternion);
		}
	}
}
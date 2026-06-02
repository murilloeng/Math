//Math
#include "Math/inc/Linear/Mat3.hpp"
#include "Math/inc/Groups/ASO3.hpp"
#include "Math/inc/Groups/ASE3.hpp"
#include "Math/inc/Groups/GSO3.hpp"
#include "Math/inc/Groups/GSE3.hpp"

namespace math
{
	namespace groups
	{
		//constructors
		GSE3::GSE3(void)
		{
			return;
		}
		GSE3::GSE3(Vec3 Vector, Quat quaternion) : m_vector{Vector}, m_quaternion{quaternion}
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
			object.m_vector_w = GSO3(m_quaternion).logarithm().Vector();
			object.m_vector_u = ASO3(object.m_vector_w).tangent_inverse().transpose() * m_vector;
			return object;
		}

		//Vector
		Vec3& GSE3::Vector(void)
		{
			return m_vector;
		}
		const Vec3& GSE3::Vector(void) const
		{
			return m_vector;
		}

		//quaternion
		Quat& GSE3::quaternion(void)
		{
			return m_quaternion;
		}
		const Quat& GSE3::quaternion(void) const
		{
			return m_quaternion;
		}

		//operators
		GSE3::operator Mat4(void) const
		{
			return Mat4(GSO3(m_quaternion), m_vector);
		}
		Vec3 GSE3::operator*(const Vec3& v) const
		{
			return GSO3(m_quaternion) * v + m_vector;
		}
		GSE3 GSE3::operator*(const GSE3& object) const
		{
			const Quat quaternion = m_quaternion * object.m_quaternion;
			const Vec3 Vector = m_vector + GSO3(m_quaternion) * object.m_vector;
			return GSE3(Vector, quaternion);
		}
	}
}
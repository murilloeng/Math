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
		ASE3::ASE3(void)
		{
			return;
		}
		ASE3::ASE3(const ASE3& object) : 
			m_vector_u(object.m_vector_u), m_vector_w(object.m_vector_w)
		{
			return;
		}
		ASE3::ASE3(const vec3& vector_u, const vec3& vector_w) : 
			m_vector_u(vector_u), m_vector_w(vector_w)
		{
			return;
		}

		//destructor
		ASE3::~ASE3(void)
		{
			return;
		}

		//matrix
		mat4 ASE3::matrix_form(void) const
		{
			return mat4(m_vector_w.spin(), m_vector_u);
		}

		//vector
		vec3& ASE3::vector_u(void)
		{
			return m_vector_u;
		}
		vec3& ASE3::vector_w(void)
		{
			return m_vector_w;
		}
		const vec3& ASE3::vector_u(void) const
		{
			return m_vector_u;
		}
		const vec3& ASE3::vector_w(void) const
		{
			return m_vector_w;
		}

		//exponential
		GSE3 ASE3::exponential(void) const
		{
			//data
			const quat q = ASO3(m_vector_w).exponential().quaternion();
			const vec3 v = ASO3(m_vector_w).tangent().transpose() * m_vector_u;
			//return
			return GSE3(v, q);
		}

		//tangent
		matrix ASE3::tangent(void) const
		{
			return matrix(6, 6);
		}
		matrix ASE3::tangent_inverse(void) const
		{
			return matrix(6, 6);
		}
		matrix ASE3::tangent_increment(void) const
		{
			return matrix(6, 6);
		}
		matrix ASE3::tangent_inverse_increment(void) const
		{
			return matrix(6, 6);
		}

		//operators
		ASE3& ASE3::operator*=(double s)
		{
			m_vector_u *= s;
			m_vector_w *= s;
			return *this;
		}
		ASE3& ASE3::operator+=(const ASE3& object)
		{
			m_vector_u += object.m_vector_u;
			m_vector_w += object.m_vector_w;
			return *this;
		}
		ASE3& ASE3::operator-=(const ASE3& object)
		{
			m_vector_u -= object.m_vector_u;
			m_vector_w -= object.m_vector_w;
			return *this;
		}

		ASE3 ASE3::operator*(double s) const
		{
			return ASE3(*this) *= s;
		}
		ASE3 ASE3::operator+(const ASE3& object) const
		{
			return ASE3(*this) += object;
		}
		ASE3 ASE3::operator-(const ASE3& object) const
		{
			return ASE3(*this) -= object;
		}

		ASE3 operator*(double s, const ASE3& object)
		{
			return ASE3(object) *= s;
		}
	}
}
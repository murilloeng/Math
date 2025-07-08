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
			//tangent
			matrix T(6, 6, mode::zeros);
			T.span(0, 0, 3, 3) = ASO3(m_vector_w).tangent();
			T.span(3, 3, 3, 3) = ASO3(m_vector_w).tangent();
			T.span(0, 3, 3, 3) = mat3(ASO3(-m_vector_w).exponential()) * ASO3(m_vector_w).tangent_increment(m_vector_u, true);
			//return
			return T;
		}
		matrix ASE3::tangent_inverse(void) const
		{
			//data
			matrix Ti(6, 6, mode::zeros);
			const mat3 Twi = ASO3(m_vector_w).tangent_inverse();
			const mat3 Awu = ASO3(m_vector_w).tangent_increment(m_vector_u, true);
			//tangent
			Ti.span(0, 0, 3, 3) = Twi;
			Ti.span(3, 3, 3, 3) = Twi;
			Ti.span(0, 3, 3, 3) = -Twi.transpose() * Awu * Twi;
			//return
			return Ti;
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
		ASE3::operator mat4(void) const
		{
			return mat4(m_vector_w.spin(), m_vector_u);
		}

		ASE3& ASE3::operator*=(double s)
		{
			m_vector_u *= s;
			m_vector_w *= s;
			return *this;
		}
		ASE3& ASE3::operator/=(double s)
		{
			m_vector_u /= s;
			m_vector_w /= s;
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
		ASE3 ASE3::operator/(double s) const
		{
			return ASE3(*this) /= s;
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
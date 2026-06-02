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
		ASE3::ASE3(void)
		{
			return;
		}
		ASE3::ASE3(const ASE3& object) : 
			m_vector_u(object.m_vector_u), m_vector_w(object.m_vector_w)
		{
			return;
		}
		ASE3::ASE3(const Vec3& vector_u, const Vec3& vector_w) : 
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
		Vec3& ASE3::vector_u(void)
		{
			return m_vector_u;
		}
		Vec3& ASE3::vector_w(void)
		{
			return m_vector_w;
		}
		const Vec3& ASE3::vector_u(void) const
		{
			return m_vector_u;
		}
		const Vec3& ASE3::vector_w(void) const
		{
			return m_vector_w;
		}

		//exponential
		GSE3 ASE3::exponential(void) const
		{
			//data
			const Quat q = ASO3(m_vector_w).exponential().quaternion();
			const Vec3 v = ASO3(m_vector_w).tangent().transpose() * m_vector_u;
			//return
			return GSE3(v, q);
		}

		//tangent
		Matrix ASE3::tangent(void) const
		{
			//tangent
			Matrix T(6, 6, mode::zeros);
			T.Span(0, 0, 3, 3) = ASO3(m_vector_w).tangent();
			T.Span(3, 3, 3, 3) = ASO3(m_vector_w).tangent();
			T.Span(0, 3, 3, 3) = Mat3(ASO3(-m_vector_w).exponential()) * ASO3(m_vector_w).tangent_increment(m_vector_u, true);
			//return
			return T;
		}
		Matrix ASE3::tangent_inverse(void) const
		{
			//data
			Matrix Ti(6, 6, mode::zeros);
			const Mat3 Twi = ASO3(m_vector_w).tangent_inverse();
			const Mat3 Awu = ASO3(m_vector_w).tangent_increment(m_vector_u, true);
			//tangent
			Ti.Span(0, 0, 3, 3) = Twi;
			Ti.Span(3, 3, 3, 3) = Twi;
			Ti.Span(0, 3, 3, 3) = -Twi.transpose() * Awu * Twi;
			//return
			return Ti;
		}
		Matrix ASE3::tangent_increment(void) const
		{
			return Matrix(6, 6);
		}
		Matrix ASE3::tangent_inverse_increment(void) const
		{
			return Matrix(6, 6);
		}

		//operators
		ASE3::operator Mat4(void) const
		{
			return Mat4(m_vector_w.spin(), m_vector_u);
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
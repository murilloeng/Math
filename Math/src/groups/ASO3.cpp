//std
#include <cmath>

//math
#include "Math/Math/inc/groups/ASO3.hpp"
#include "Math/Math/inc/groups/GSO3.hpp"

namespace math
{
	namespace groups
	{
		//constructor
		ASO3::ASO3(void)
		{
			return;
		}
		ASO3::ASO3(vec3 vector) : m_vector{vector}
		{
			return;
		}

		//destructor
		ASO3::~ASO3(void)
		{
			return;
		}

		//matrix
		mat3 ASO3::matrix(void) const
		{
			//data
			mat3 matrix;
			const double* p = m_vector.data();
			//matrix
			matrix[7] = -p[0];
			matrix[2] = -p[1];
			matrix[3] = -p[2];
			matrix[5] = +p[0];
			matrix[6] = +p[1];
			matrix[1] = +p[2];
			matrix[0] = matrix[4] = matrix[8] = 0;
			//return
			return matrix;
		}

		//vector
		vec3& ASO3::vector(void)
		{
			return m_vector;
		}
		const vec3& ASO3::vector(void) const
		{
			return m_vector;
		}

		//exponential
		GSO3 ASO3::exponential(void) const
		{
			//group
			GSO3 group;
			const double t = m_vector.norm();
			group.m_quaternion[0] = cos(t / 2);
			group.m_quaternion[1] = t ? sin(t / 2) * m_vector[0] / t : 0;
			group.m_quaternion[2] = t ? sin(t / 2) * m_vector[1] / t : 0;
			group.m_quaternion[3] = t ? sin(t / 2) * m_vector[2] / t : 0;
			//return
			return group;
		}

		//tangent
		mat3 ASO3::tangent(void) const
		{
			const mat3 s = m_vector.spin();
			const double t = m_vector.norm();
			return mat3::eye() - fn(t, 2) * s + fn(t, 3) * s * s;
		}
		vec3 ASO3::tangent(const vec3& a) const
		{
			const vec3 b = m_vector.cross(a);
			const vec3 c = m_vector.cross(b);
			const double t = m_vector.norm();
			return a - fn(t, 2) * b + fn(t, 3) * c;
		}

		mat3 ASO3::tangent_inverse(void) const
		{
			const mat3 s = m_vector.spin();
			const double t = m_vector.norm();
			return mat3::eye() + s / 2 + (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2 * s * s;
		}
		vec3 ASO3::tangent_inverse(const vec3& a) const
		{
			const vec3 b = m_vector.cross(a);
			const vec3 c = m_vector.cross(b);
			const double t = m_vector.norm();
			return a + b / 2 + (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2 * c;
		}

		mat3 ASO3::tangent_increment(const vec3& a) const
		{
			const vec3& x = m_vector;
			const vec3 b = m_vector.cross(a);
			const vec3 c = m_vector.cross(b);
			const double t = m_vector.norm();
			return -dfn(t, 2) / t * b.outer(x) + fn(t, 2) * a.spin() + 
			dfn(t, 3) / t * c.outer(x) - fn(t, 3) * (x.spin() * a.spin() + b.spin());
		}
		vec3 ASO3::tangent_increment(const vec3& a, const vec3& u) const
		{
			const vec3& x = m_vector;
			const vec3 b = m_vector.cross(a);
			const vec3 c = m_vector.cross(b);
			const double t = m_vector.norm();
			const double s = m_vector.inner(u);
			return -dfn(t, 2) / t * s * b + fn(t, 2) * a.cross(u) + 
			dfn(t, 3) / t * s * c - fn(t, 3) * (x.spin() * a.cross(u) + b.cross(u));
		}

		mat3 ASO3::tangent_inverse_increment(const vec3& a) const
		{
			const mat3 s = m_vector.spin();
			const double t = m_vector.norm();
			return mat3::eye() + s / 2 + (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2 * s * s;
		}
		vec3 ASO3::tangent_inverse_increment(const vec3& a, const vec3& u) const
		{
			return vec3();
		}

		//operators
		ASO3& ASO3::operator*=(double s)
		{
			m_vector *= s;
			return *this;
		}
		ASO3& ASO3::operator+=(const ASO3& a)
		{
			m_vector += a.m_vector;
			return *this;
		}
		ASO3& ASO3::operator-=(const ASO3& a)
		{
			m_vector -= a.m_vector;
			return *this;
		}

		ASO3 ASO3::operator*(double s) const
		{
			return s * m_vector;
		}
		ASO3 ASO3::operator+(const ASO3& a) const
		{
			return m_vector + a.m_vector;
		}
		ASO3 ASO3::operator-(const ASO3& a) const
		{
			return m_vector - a.m_vector;
		}

		//rotation
		double ASO3::fn(double t, uint32_t n)
		{
			if(fabs(t) < 2 * M_PI)
			{
				return funt(t, n);
			}
			else
			{
				if(n % 2)
				{
					const uint32_t p = n / 2;
					return (p % 2 ? -1 : +1) / pow(t, n) * (sin(t) - sint(t, p));
				}
				else
				{
					const uint32_t p = n / 2;
					return (p % 2 ? -1 : +1) / pow(t, n) * (cos(t) - cost(t, p));
				}
			}
		}
		double ASO3::dfn(double t, uint32_t n)
		{
			return t * (n * fn(t, n + 2) - fn(t, n + 1));
		}
		double ASO3::funt(double t, uint32_t n)
		{
			//data
			int32_t s = 1;
			uint32_t k = 0;
			double v = 0, z = 1, dv;
			double p = std::tgamma(n + 1);
			//compute
			while(true)
			{
				dv = s * z / p;
				if(v + dv == v)
				{
					return v;
				}
				else
				{
					k++;
					v += dv;
					s *= -1;
					z *= t * t;
					p *= 2 * k + n;
					p *= 2 * k + n - 1;
				}
			}
		}
		double ASO3::cost(double t, uint32_t n)
		{
			int32_t s = 1;
			uint32_t a = 1;
			double v = 0, p = 1;
			for(uint32_t k = 0; k < n; k++)
			{
				v += s * p / a;
				s *= -1;
				p *= t * t;
				a *= 2 * k + 1;
				a *= 2 * k + 2;
			}
			return v;
		}
		double ASO3::sint(double t, uint32_t n)
		{
			int32_t s = 1;
			uint32_t a = 1;
			double v = 0, p = t;
			for(uint32_t k = 0; k < n; k++)
			{
				v += s * p / a;
				s *= -1;
				p *= t * t;
				a *= 2 * k + 2;
				a *= 2 * k + 3;
			}
			return v;
		}
	}
}
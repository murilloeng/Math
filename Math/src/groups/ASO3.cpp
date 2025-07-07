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

		}

		//destructor
		ASO3::~ASO3(void)
		{
			return;
		}

		//exponential
		GSO3 ASO3::exponential(void) const
		{
			return GSO3();
		}

		//tangent
		mat3 ASO3::tangent(void) const
		{
			return mat3();
		}
		mat3 ASO3::tangent_inverse(void) const
		{
			return mat3();
		}
		mat3 ASO3::tangent_increment(const vec3&) const
		{
			return mat3();
		}
		mat3 ASO3::tangent_inverse_increment(const vec3&) const
		{
			return mat3();
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
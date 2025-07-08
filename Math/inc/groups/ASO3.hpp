#pragma once

//math
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat3.hpp"

namespace math
{
	namespace groups
	{
		class GSO3;
	}
}

namespace math
{
	namespace groups
	{
		class ASO3
		{
		public:
			//constructor
			ASO3(void);
			ASO3(vec3);

			//destructor
			~ASO3(void);

			//vector
			vec3& vector(void);
			const vec3& vector(void) const;

			//exponential
			GSO3 exponential(void) const;

			//tangent
			mat3 tangent(bool = false) const;
			vec3 tangent(const vec3&, bool = false) const;

			mat3 tangent_inverse(bool = false) const;
			vec3 tangent_inverse(const vec3&, bool = false) const;

			mat3 tangent_increment(const vec3&, bool = false) const;
			vec3 tangent_increment(const vec3&, const vec3&, bool = false) const;

			mat3 tangent_inverse_increment(const vec3&, bool = false) const;
			vec3 tangent_inverse_increment(const vec3&, const vec3&, bool = false) const;

			//operators
			operator mat3(void) const;

			ASO3& operator*=(double);
			ASO3& operator/=(double);
			ASO3& operator+=(const ASO3&);
			ASO3& operator-=(const ASO3&);

			ASO3 operator*(double) const;
			ASO3 operator/(double) const;
			ASO3 operator+(const ASO3&) const;
			ASO3 operator-(const ASO3&) const;

			friend ASO3 operator*(double, const ASO3&);

		private:
			//rotation
			static double fn(double, uint32_t);
			static double dfn(double, uint32_t);
			static double funt(double, uint32_t);
			static double cost(double, uint32_t);
			static double sint(double, uint32_t);

			//data
			vec3 m_vector;
		};
	}
}
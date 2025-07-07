#pragma once

//math
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat4.hpp"

namespace math
{
	namespace groups
	{
		class GSE3;
	}
}

namespace math
{
	namespace groups
	{
		class ASE3
		{
		public:
			//constructors
			ASE3(void);
			ASE3(const ASE3&);
			ASE3(const vec3&, const vec3&);

			//destructor
			~ASE3(void);

			//matrix
			mat4 matrix_form(void) const;

			//vector
			vec3& vector_u(void);
			vec3& vector_w(void);
			const vec3& vector_u(void) const;
			const vec3& vector_w(void) const;

			//exponential
			GSE3 exponential(void) const;

			//tangent
			matrix tangent(void) const;
			matrix tangent_inverse(void) const;
			matrix tangent_increment(void) const;
			matrix tangent_inverse_increment(void) const;

			//operators
			ASE3& operator*=(double);
			ASE3& operator+=(const ASE3&);
			ASE3& operator-=(const ASE3&);

			ASE3 operator*(double) const;
			ASE3 operator+(const ASE3&) const;
			ASE3 operator-(const ASE3&) const;

			friend ASE3 operator*(double, const ASE3&);

		private:
			//data
			vec3 m_vector_u;
			vec3 m_vector_w;

			//friends
			friend class GSE3;
		};
	}
}
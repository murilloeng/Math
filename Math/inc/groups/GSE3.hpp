#pragma once

//math
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/quat.hpp"
#include "Math/Math/inc/linear/mat4.hpp"

namespace math
{
	namespace groups
	{
		class ASE3;
	}
}

namespace math
{
	namespace groups
	{
		class GSE3
		{
		public:
			//constructors
			GSE3(void);
			GSE3(vec3, quat);

			//destructor
			~GSE3(void);
			
			//inverse
			GSE3 inverse(void) const;
			
			//logarithm
			ASE3 logarithm(void) const;

			//matrix
			mat4 matrix_form(void) const;

			//vector
			vec3& vector(void);
			const vec3& vector(void) const;

			//quaternion
			quat& quaternion(void);
			const quat& quaternion(void) const;

			//operators
			vec3 operator*(const vec3&) const;
			GSE3 operator*(const GSE3&) const;

		private:
			//data
			vec3 m_vector;
			quat m_quaternion;

			//friends
			friend class ASE3;
		};
	}
}
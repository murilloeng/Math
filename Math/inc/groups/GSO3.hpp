#pragma once

//math
#include "Math/Math/inc/linear/quat.hpp"

namespace math
{
	namespace groups
	{
		class ASO3;
	}
}

namespace math
{
	namespace groups
	{
		class GSO3
		{
		public:
			//constructors
			GSO3(void);
			GSO3(quat);

			//destructor
			~GSO3(void);
			
			//inverse
			GSO3 inverse(void) const;
			
			//logarithm
			ASO3 logarithm(void) const;

			//matrix
			mat3 matrix_form(void) const;

			//quaternion
			quat& quaternion(void);
			const quat& quaternion(void) const;

			//operators
			vec3 operator*(const vec3&) const;
			GSO3 operator*(const GSO3&) const;

		private:
			//data
			quat m_quaternion;

			//friends
			friend class ASO3;
		};
	}
}
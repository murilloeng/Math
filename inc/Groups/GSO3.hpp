#pragma once

//Math
#include "Math/inc/Linear/Quat.hpp"

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
			GSO3(Quat);

			//destructor
			~GSO3(void);

			//inverse
			GSO3 inverse(void) const;

			//logarithm
			ASO3 logarithm(void) const;

			//quaternion
			Quat& quaternion(void);
			const Quat& quaternion(void) const;

			//operators
			operator Mat3(void) const;
			Vec3 operator*(const Vec3&) const;
			GSO3 operator*(const GSO3&) const;

		private:
			//data
			Quat m_quaternion;

			//friends
			friend class ASO3;
		};
	}
}
#pragma once

//Math
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Quat.hpp"
#include "Math/inc/Linear/Mat4.hpp"

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
			GSE3(Vec3, Quat);

			//destructor
			~GSE3(void);

			//inverse
			GSE3 inverse(void) const;

			//logarithm
			ASE3 logarithm(void) const;

			//vector
			Vec3& Vector(void);
			const Vec3& Vector(void) const;

			//quaternion
			Quat& quaternion(void);
			const Quat& quaternion(void) const;

			//operators
			operator Mat4(void) const;
			Vec3 operator*(const Vec3&) const;
			GSE3 operator*(const GSE3&) const;

		private:
			//data
			Vec3 m_vector;
			Quat m_quaternion;

			//friends
			friend class ASE3;
		};
	}
}
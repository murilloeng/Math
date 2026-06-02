#pragma once

//Math
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Mat3.hpp"

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
			ASO3(Vec3);

			//destructor
			~ASO3(void);

			//Vector
			Vec3& Vector(void);
			const Vec3& Vector(void) const;

			//exponential
			GSO3 exponential(void) const;

			//tangent
			Mat3 tangent(bool = false) const;
			Vec3 tangent(const Vec3&, bool = false) const;

			Mat3 tangent_inverse(bool = false) const;
			Vec3 tangent_inverse(const Vec3&, bool = false) const;

			Mat3 tangent_increment(const Vec3&, bool = false) const;
			Vec3 tangent_increment(const Vec3&, const Vec3&, bool = false) const;

			Mat3 tangent_inverse_increment(const Vec3&, bool = false) const;
			Vec3 tangent_inverse_increment(const Vec3&, const Vec3&, bool = false) const;

			//operators
			operator Mat3(void) const;

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
			Vec3 m_vector;
		};
	}
}
#pragma once

//Math
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Mat4.hpp"

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
			ASE3(const Vec3&, const Vec3&);

			//destructor
			~ASE3(void);

			//Vector
			Vec3& vector_u(void);
			Vec3& vector_w(void);
			const Vec3& vector_u(void) const;
			const Vec3& vector_w(void) const;

			//exponential
			GSE3 exponential(void) const;

			//tangent
			Matrix tangent(void) const;
			Matrix tangent_inverse(void) const;
			Matrix tangent_increment(void) const;
			Matrix tangent_inverse_increment(void) const;

			//operators
			operator Mat4(void) const;

			ASE3& operator*=(double);
			ASE3& operator/=(double);
			ASE3& operator+=(const ASE3&);
			ASE3& operator-=(const ASE3&);

			ASE3 operator*(double) const;
			ASE3 operator/(double) const;
			ASE3 operator+(const ASE3&) const;
			ASE3 operator-(const ASE3&) const;

			friend ASE3 operator*(double, const ASE3&);

		private:
			//data
			Vec3 m_vector_u;
			Vec3 m_vector_w;

			//friends
			friend class GSE3;
		};
	}
}
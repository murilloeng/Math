#pragma once

//std
#include <cstdint>

namespace math
{
	namespace linear
	{
		class Mat3;
	}
}

namespace math
{
	namespace linear
	{
		class Vec3
		{
		public:
			//constructors
			Vec3(const Vec3&) = default;
			explicit Vec3(const double*);
			explicit Vec3(double = 0, double = 0, double = 0);

			//destructor
			~Vec3(void) = default;

			//data
			double* data(void);
			const double* data(void) const;

			//operators
			Vec3 operator+(void) const;
			Vec3 operator-(void) const;

			Vec3 operator/(double) const;
			Vec3 operator+(const Vec3&) const;
			Vec3 operator-(const Vec3&) const;

			Vec3& operator+=(double);
			Vec3& operator-=(double);
			Vec3& operator*=(double);
			Vec3& operator/=(double);

			Vec3& operator+=(const Vec3&);
			Vec3& operator-=(const Vec3&);

			Vec3& operator=(const Vec3&) = default;

			double& operator[](uint32_t);
			const double& operator[](uint32_t) const;

			//norm
			double norm(void) const;

			//direction
			Vec3 unit(void) const;
			Vec3& normalize(void);
			Vec3& project(const Vec3&);

			//projection
			Mat3 spin(void) const;
			Mat3 projection(void) const;

			//products
			Mat3 outer(const Vec3&) const;
			Vec3 cross(const Vec3&) const;
			double inner(const Vec3&) const;

			//interpolation
			Vec3 lerp(const Vec3&, double) const;

			//friends
			friend Vec3 operator*(double, const Vec3&);

		private:
			//data
			double m_data[3];

			//friends
			friend class math::linear::Mat3;
		};
	}
}
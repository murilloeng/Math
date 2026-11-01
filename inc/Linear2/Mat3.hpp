#pragma once

//std
#include <cstdint>
#include <initializer_list>

namespace math
{
	namespace linear
	{
		class Vec3;
	}
}

namespace math
{
	namespace linear
	{
		class Mat3
		{
		public:
			//constructors
			Mat3(void);
			Mat3(const Mat3&) = default;
			explicit Mat3(const double*);
			explicit Mat3(const Vec3&, const Vec3&, const Vec3&);
			Mat3(std::initializer_list<std::initializer_list<double>>);

			//destructor
			~Mat3(void) = default;

			//data
			double* data(void);
			const double* data(void) const;

			//operators
			Mat3 operator+(void) const;
			Mat3 operator-(void) const;

			Mat3 operator/(double) const;
			Mat3 operator+(const Mat3&) const;
			Mat3 operator-(const Mat3&) const;
			Mat3 operator*(const Mat3&) const;
			Vec3 operator*(const Vec3&) const;

			Mat3& operator+=(double);
			Mat3& operator-=(double);
			Mat3& operator*=(double);
			Mat3& operator/=(double);

			Mat3& operator+=(const Mat3&);
			Mat3& operator-=(const Mat3&);

			double& operator[](uint32_t);
			double& operator()(uint32_t, uint32_t);

			const double& operator[](uint32_t) const;
			const double& operator()(uint32_t, uint32_t) const;

			//linear
			static Mat3 eye(void);
			Mat3 inverse(void) const;
			Mat3 transpose(void) const;
			Vec3 eigenvalues(void) const;
			Mat3 eigenvectors(void) const;

			double trace(void) const;
			double lode_angle(void) const;
			double determinant(void) const;
			double invariant_I1(void) const;
			double invariant_I2(void) const;
			double invariant_I3(void) const;
			double deviatoric_J2(void) const;
			double deviatoric_J3(void) const;

			//friends
			friend Mat3 operator*(double, const Mat3&);

		private:
			//data
			double m_data[9];

			//friends
			friend class math::linear::Vec3;
		};
	}
}
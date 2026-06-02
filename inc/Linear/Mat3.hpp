#pragma once

//Math
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	class Vec3;
	class Quat;
}

namespace math
{
	class Mat3 : public Matrix
	{
	public:
		//constructors
		Mat3(double*);
		Mat3(const Mat3&);
		Mat3(const double*);
		Mat3(mode = mode::null);
		Mat3(std::initializer_list<double>);
		Mat3(const Vec3&, const Vec3&, const Vec3&);
		Mat3(const double*, const double*, const double*);

		//destructor
		~Mat3(void);

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

		Mat3& operator=(const Mat3&);

		Mat3& operator+=(const Mat3&);
		Mat3& operator-=(const Mat3&);

		//linear
		static Mat3 eye(void);
		Vec3 eigen(void) const;
		Mat3 inverse(void) const;
		Mat3 transpose(void) const;

		double trace(void) const;
		double lode_angle(void) const;
		double determinant(void) const;
		double invariant_I1(void) const;
		double invariant_I2(void) const;
		double invariant_I3(void) const;
		double deviatoric_J2(void) const;
		double deviatoric_J3(void) const;

		//rotation
		Vec3 rotation(void) const;
		Quat quaternion(void) const;

		//friends
		friend Mat3 operator*(double, const Mat3&);
	};
}
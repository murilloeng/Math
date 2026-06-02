#pragma once

//Math
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	class Vec3;
}

namespace math
{
	class Mat4 : public Matrix
	{
	public:
		//constructors
		Mat4(double*);
		Mat4(const Mat4&);
		Mat4(const double*);
		Mat4(mode = mode::null);
		Mat4(const Mat3&, const Vec3&);
		Mat4(std::initializer_list<double>);

		//destructor
		~Mat4(void);

		//operators
		Mat4 operator+(void) const;
		Mat4 operator-(void) const;

		Mat4 operator/(double) const;
		Mat4 operator+(const Mat4&) const;
		Mat4 operator-(const Mat4&) const;
		Mat4 operator*(const Mat4&) const;
		Vec3 operator*(const Vec3&) const;

		Mat4& operator+=(double);
		Mat4& operator-=(double);
		Mat4& operator*=(double);
		Mat4& operator/=(double);

		Mat4& operator=(const Mat4&);

		Mat4& operator+=(const Mat4&);
		Mat4& operator-=(const Mat4&);

		//linear
		static Mat4 eye(void);
		Mat4 transpose(void) const;

		//friends
		friend Mat4 operator*(double, const Mat4&);
	};
}
#pragma once

//std
#include <cstdint>

//Math
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	class Vec2;
}

namespace math
{
	class Mat2 : public Matrix
	{
	public:
		//constructors
		Mat2(double*);
		Mat2(const Mat2&);
		Mat2(const double*);
		Mat2(mode = mode::null);
		Mat2(const Vec2&, const Vec2&);
		Mat2(const double*, const double*);
		Mat2(std::initializer_list<double>);

		//destructor
		~Mat2(void);

		//operators
		Mat2 operator+(void) const;
		Mat2 operator-(void) const;

		Mat2 operator/(double) const;
		Mat2 operator+(const Mat2&) const;
		Mat2 operator-(const Mat2&) const;
		Mat2 operator*(const Mat2&) const;
		Vec2 operator*(const Vec2&) const;

		Mat2& operator+=(double);
		Mat2& operator-=(double);
		Mat2& operator*=(double);
		Mat2& operator/=(double);

		Mat2& operator=(const Mat2&);

		Mat2& operator+=(const Mat2&);
		Mat2& operator-=(const Mat2&);

		double& operator[](uint32_t);
		double& operator()(uint32_t);
		double& operator()(uint32_t, uint32_t);

		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t) const;
		const double& operator()(uint32_t, uint32_t) const;

		//linear
		static Mat2 eye(void);
		Vec2 eigen(void) const;
		Mat2 inverse(void) const;
		Mat2 transpose(void) const;

		double trace(void) const;
		double determinant(void) const;
		double invariant_I1(void) const;
		double invariant_I2(void) const;
		double deviatoric_J2(void) const;

		//friends
		friend Mat2 operator*(double, const Mat2&);
	};
}
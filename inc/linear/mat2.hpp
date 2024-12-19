#pragma once

//std
#include <cstdint>

//math
#include "Math/inc/linear/matrix.hpp"

namespace math
{
	class vec2;
}

namespace math
{
	class mat2 : public matrix
	{
	public:
		//constructors
		mat2(double*);
		mat2(const mat2&);
		mat2(const double*);
		mat2(mode = mode::null);
		mat2(const vec2&, const vec2&);
		mat2(const double*, const double*);
		mat2(std::initializer_list<double>);

		//destructor
		virtual ~mat2(void);

		//operators
		mat2 operator+(void) const;
		mat2 operator-(void) const;

		mat2 operator/(double) const;
		mat2 operator+(const mat2&) const;
		mat2 operator-(const mat2&) const;
		mat2 operator*(const mat2&) const;
		vec2 operator*(const vec2&) const;

		mat2& operator+=(double);
		mat2& operator-=(double);
		mat2& operator*=(double);
		mat2& operator/=(double);

		mat2& operator=(const mat2&);

		mat2& operator+=(const mat2&);
		mat2& operator-=(const mat2&);

		double& operator[](uint32_t);
		double& operator()(uint32_t);
		double& operator()(uint32_t, uint32_t);

		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t) const;
		const double& operator()(uint32_t, uint32_t) const;

		//linear
		static mat2 eye(void);
		vec2 eigen(void) const;
		mat2 inverse(void) const;
		mat2 transpose(void) const;

		double trace(void) const;
		double determinant(void) const;
		double invariant_I1(void) const;
		double invariant_I2(void) const;
		double deviatoric_J2(void) const;

		//friends
		friend mat2 operator*(double, const mat2&);
	};
}
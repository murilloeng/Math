#pragma once

//math
#include "Math/Math/inc/linear/matrix.hpp"

namespace math
{
	class vec3;
}

namespace math
{
	class mat4 : public matrix
	{
	public:
		//constructors
		mat4(double*);
		mat4(const mat4&);
		mat4(const double*);
		mat4(mode = mode::null);
		mat4(std::initializer_list<double>);

		//destructor
		~mat4(void);

		//operators
		mat4 operator+(void) const;
		mat4 operator-(void) const;

		mat4 operator/(double) const;
		mat4 operator+(const mat4&) const;
		mat4 operator-(const mat4&) const;
		mat4 operator*(const mat4&) const;
		vec3 operator*(const vec3&) const;

		mat4& operator+=(double);
		mat4& operator-=(double);
		mat4& operator*=(double);
		mat4& operator/=(double);

		mat4& operator=(const mat4&);

		mat4& operator+=(const mat4&);
		mat4& operator-=(const mat4&);

		//linear
		static mat4 eye(void);
		mat4 transpose(void) const;

		//friends
		friend mat4 operator*(double, const mat4&);
	};
}
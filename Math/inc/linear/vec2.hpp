#pragma once

//math
#include "Galileo/mat/inc/linear/vector.hpp"

namespace math
{
	class vec2 : public vector
	{
	public:
		//constructors
		vec2(void);
		vec2(double*);
		vec2(const vec2&);
		vec2(const double*);
		vec2(double, double);

		//destructor
		~vec2(void);

		//operators
		vec2 operator+(void) const;
		vec2 operator-(void) const;

		vec2 operator/(double) const;
		vec2 operator+(const vec2&) const;
		vec2 operator-(const vec2&) const;

		vec2& operator+=(double);
		vec2& operator-=(double);
		vec2& operator*=(double);
		vec2& operator/=(double);

		vec2& operator=(const vec2&);

		vec2& operator+=(const vec2&);
		vec2& operator-=(const vec2&);

		double& operator[](uint32_t);
		double& operator()(uint32_t);

		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t) const;

		//linear
		vec2& normalize(void);
		vec2& project(const vec2&);

		double inner(const vec2&) const;
		double cross(const vec2&) const;

		//friends
		friend vec2 operator*(double, const vec2&);
	};
}
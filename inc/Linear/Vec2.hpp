#pragma once

//Math
#include "Math/inc/Linear/Vector.hpp"

namespace math
{
	class Vec2 : public Vector
	{
	public:
		//constructors
		Vec2(void);
		Vec2(double*);
		Vec2(const Vec2&);
		Vec2(const double*);
		Vec2(double, double);

		//destructor
		~Vec2(void);

		//operators
		Vec2 operator+(void) const;
		Vec2 operator-(void) const;

		Vec2 operator/(double) const;
		Vec2 operator+(const Vec2&) const;
		Vec2 operator-(const Vec2&) const;

		Vec2& operator+=(double);
		Vec2& operator-=(double);
		Vec2& operator*=(double);
		Vec2& operator/=(double);

		Vec2& operator=(const Vec2&);

		Vec2& operator+=(const Vec2&);
		Vec2& operator-=(const Vec2&);

		double& operator[](uint32_t);
		double& operator()(uint32_t);

		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t) const;

		//linear
		Vec2& normalize(void);
		Vec2& project(const Vec2&);

		double inner(const Vec2&) const;
		double cross(const Vec2&) const;

		//friends
		friend Vec2 operator*(double, const Vec2&);
	};
}
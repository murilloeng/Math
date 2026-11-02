#pragma once

//std
#include <cstdint>
#include <initializer_list>

namespace math
{
	class Mat2;
}

namespace math
{
	class Vec2
	{
	public:
		//constructors
		Vec2(const double*);
		Vec2(const Vec2&) = default;
		Vec2(double = 0, double = 0);
		Vec2(std::initializer_list<double>);

		//destructor
		~Vec2(void) = default;

		//data
		double* data(void);
		const double* data(void) const;

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

		Vec2& operator+=(const Vec2&);
		Vec2& operator-=(const Vec2&);

		Vec2& operator=(const Vec2&) = default;

		double& operator[](uint32_t);
		const double& operator[](uint32_t) const;

		//norm
		double norm(void) const;

		//direction
		Vec2 unit(void) const;
		Vec2& normalize(void);

		//projection
		Mat2 projection(void) const;

		//products
		double inner(const Vec2&) const;
		double cross(const Vec2&) const;

		//interpolation
		Vec2 lerp(const Vec2&, double) const;

		//friends
		friend Vec2 operator*(double, const Vec2&);

	private:
		//data
		double m_data[2];

		//friends
		friend class math::Mat2;
	};
}
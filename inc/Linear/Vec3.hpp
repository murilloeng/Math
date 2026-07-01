#pragma once

//Math
#include "Math/inc/Linear/Vector.hpp"

namespace math
{
	class Quat;
	class Mat3;
}

namespace math
{
	class Vec3 : public Vector
	{
	public:
		//constructors
		Vec3(void);
		Vec3(double*);
		Vec3(const Vec3&);
		Vec3(const double*);
		Vec3(double, double, double);

		//destructor
		~Vec3(void);

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

		Vec3& operator=(const Vec3&);

		Vec3& operator+=(const Vec3&);
		Vec3& operator-=(const Vec3&);
		Vec3& operator*=(const Mat3&);

		//linear
		Vec3 unit(void) const;
		Vec3& normalize(void);
		Vec3& project(const Vec3&);
		const Vec3& triad(Vec3&, Vec3&, double = 0) const;

		Mat3 spin(void) const;
		Quat quaternion(void) const;
		Mat3 projection(void) const;

		Mat3 outer(const Vec3&) const;
		Vec3 cross(const Vec3&) const;
		Vec3 rotate(const Vec3&) const;
		double inner(const Vec3&) const;

		Vec3 lerp(const Vec3&, double) const;

		//static
		static Vec3 base(uint32_t);

		//rotation
		Mat3 rotation_tensor(void) const;

		Mat3 rotation_gradient(bool = false) const;
		Vec3 rotation_gradient(const Vec3&, bool = false) const;

		Mat3 rotation_class(uint32_t, bool = false) const;
		Vec3 rotation_class(const Vec3&, uint32_t, bool = false) const;

		Mat3 rotation_gradient_inverse(bool = false) const;
		Vec3 rotation_gradient_inverse(const Vec3&, bool = false) const;

		Mat3 rotation_hessian(uint32_t, bool = false) const;
		Mat3 rotation_hessian(const Vec3&, bool = false) const;
		Vec3 rotation_hessian(const Vec3&, const Vec3&, bool = false) const;

		Mat3 rotation_hessian_inverse(uint32_t, bool = false) const;
		Mat3 rotation_hessian_inverse(const Vec3&, bool = false) const;
		Vec3 rotation_hessian_inverse(const Vec3&, const Vec3&, bool = false) const;

		Mat3 rotation_class_increment(const Vec3&, uint32_t, bool = false) const;
		Vec3 rotation_class_increment(const Vec3&, const Vec3&, uint32_t, bool = false) const;

		Mat3 rotation_third(uint32_t, uint32_t, bool = false) const;
		Mat3 rotation_third(const Vec3&, const Vec3&, bool = false, bool = true) const;
		Mat3 rotation_third_inverse(const Vec3&, const Vec3&, bool = false, bool = true) const;

		//friends
		friend Vec3 operator*(double, const Vec3&);
	};
}
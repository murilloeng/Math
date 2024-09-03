#pragma once

//math
#include "Math/inc/linear/vector.hpp"

namespace math
{
	class quat;
	class mat3;
}

namespace math
{
	class vec3 : public vector
	{
	public:
		//constructors
		vec3(void);
		vec3(double*);
		vec3(const vec3&);
		vec3(const double*);
		vec3(double, double, double);

		//destructor
		virtual ~vec3(void);

		//operators
		vec3 operator+(void) const;
		vec3 operator-(void) const;

		vec3 operator/(double) const;
		vec3 operator+(const vec3&) const;
		vec3 operator-(const vec3&) const;

		vec3& operator+=(double);
		vec3& operator-=(double);
		vec3& operator*=(double);
		vec3& operator/=(double);

		vec3& operator=(const vec3&);

		vec3& operator+=(const vec3&);
		vec3& operator-=(const vec3&);
		vec3& operator*=(const mat3&);

		double& operator[](uint32_t);
		double& operator()(uint32_t);

		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t) const;

		//linear
		vec3& normalize(void);
		vec3& project(const vec3&);
		const vec3& triad(vec3&, vec3&, double = 0) const;

		mat3 spin(void) const;
		quat quaternion(void) const;
		mat3 projection(void) const;

		mat3 outer(const vec3&) const;
		vec3 cross(const vec3&) const;
		vec3 rotate(const vec3&) const;
		double inner(const vec3&) const;

		vec3 lerp(const vec3&, double) const;

		mat3 rotation_tensor(void) const;

		mat3 rotation_gradient(bool = false) const;
		vec3 rotation_gradient(const vec3&, bool = false) const;

		mat3 rotation_class(uint32_t, bool = false) const;
		vec3 rotation_class(const vec3&, uint32_t, bool = false) const;

		mat3 rotation_gradient_inverse(bool = false) const;
		vec3 rotation_gradient_inverse(const vec3&, bool = false) const;

		mat3 rotation_hessian(const vec3&, bool = false) const;
		vec3 rotation_hessian(const vec3&, const vec3&, bool = false) const;

		mat3 rotation_hessian_inverse(const vec3&, bool = false) const;
		vec3 rotation_hessian_inverse(const vec3&, const vec3&, bool = false) const;

		mat3 rotation_class_increment(const vec3&, uint32_t, bool = false) const;
		vec3 rotation_class_increment(const vec3&, const vec3&, uint32_t, bool = false) const;

		mat3 rotation_higher(const vec3&, const vec3&, bool = false, bool = true) const;
		mat3 rotation_higher_inverse(const vec3&, const vec3&, bool = false, bool = true) const;

		//friends
		friend vec3 operator*(double, const vec3&);
	};
}
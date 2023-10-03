#pragma once

//mat
#include "inc/linear/matrix.hpp"

namespace math
{
	class vec3;
	class quat;
}

namespace math
{
	class mat3 : public matrix
	{
	public:
		//constructors
		mat3(void);
		mat3(double*);
		mat3(const mat3&);
		mat3(const double*);
		mat3(std::initializer_list<double>);
		mat3(const vec3&, const vec3&, const vec3&);
		mat3(const double*, const double*, const double*);

		//destructor
		virtual ~mat3(void);

		//operators
		mat3 operator+(void) const;
		mat3 operator-(void) const;

		mat3 operator/(double) const;
		mat3 operator+(const mat3&) const;
		mat3 operator-(const mat3&) const;
		mat3 operator*(const mat3&) const;
		vec3 operator*(const vec3&) const;

		mat3& operator+=(double);
		mat3& operator-=(double);
		mat3& operator*=(double);
		mat3& operator/=(double);

		mat3& operator=(const mat3&);

		mat3& operator+=(const mat3&);
		mat3& operator-=(const mat3&);

		double& operator[](unsigned);
		double& operator()(unsigned);
		double& operator()(unsigned, unsigned);

		const double& operator[](unsigned) const;
		const double& operator()(unsigned) const;
		const double& operator()(unsigned, unsigned) const;

		//misc
		static mat3 eye(void);

		//linear
		vec3 eigen(void) const;
		mat3 inverse(void) const;
		mat3 transpose(void) const;

		double trace(void) const;
		double lode_angle(void) const;
		double determinant(void) const;
		double invariant_I1(void) const;
		double invariant_I2(void) const;
		double invariant_I3(void) const;
		double deviatoric_J2(void) const;
		double deviatoric_J3(void) const;

		//rotation
		vec3 rotation(void) const;
		quat quaternion(void) const;

		//friends
		friend mat3 operator*(double, const mat3&);
	};
}
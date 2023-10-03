#pragma once

//math
#include "inc/linear/vector.hpp"

namespace math
{
	class vec3;
	class mat3;
}

namespace math
{
	class quat : public vector
	{
	public:
		//constructors
		quat(void);
		quat(double*);
		quat(const quat&);
		quat(const double*);
		quat(double, unsigned);
		quat(double, const vec3&);
		quat(double, double, double, double);
		quat(const vec3&, const vec3&, const vec3&);
		quat(const double*, const double*, const double*);

		//destructor
		virtual ~quat(void);

		//operators
		quat& operator=(const quat&);
		quat& operator*=(const quat&);

		quat operator*(const quat&) const;

		double& operator[](unsigned);
		double& operator()(unsigned);

		const double& operator[](unsigned) const;
		const double& operator()(unsigned) const;

		//views
		quat& reset(void);
		quat& view_x(void);
		quat& view_y(void);
		quat& view_z(void);
		quat& isometric(unsigned);

		//linear
		vec3 axial(void) const;
		vec3 pseudo(void) const;
		vec3 pseudo(vec3) const;
		double angle(void) const;

		quat conjugate(void) const;

		vec3 rotate(const vec3&) const;
		vec3 conjugate(const vec3&) const;
		quat conjugate(const quat&) const;

		quat slerp(const quat&, double) const;

		mat3 rotation(void) const;
		double* rotation(double*) const;
	};
}
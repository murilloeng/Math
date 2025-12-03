#pragma once

//std
#include <cstdint>

//math
#include "Math/Math/inc/linear/vector.hpp"

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
		quat(double, uint32_t);
		quat(double, const vec3&);
		quat(double, double, double, double);
		quat(const vec3&, const vec3&, const vec3&);
		quat(const double*, const double*, const double*);

		//destructor
		~quat(void);

		//operators
		quat& operator=(const quat&);
		quat& operator*=(const quat&);

		quat operator*(const quat&) const;

		double& operator[](uint32_t);
		double& operator()(uint32_t);

		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t) const;

		//views
		quat& reset(void);
		quat& view_x(void);
		quat& view_y(void);
		quat& view_z(void);
		quat& isometric(uint32_t);

		//linear
		quat& normalize(void);
		vec3 axial(void) const;
		vec3 pseudo(void) const;
		double angle(void) const;
		vec3 pseudo(const vec3&) const;
		vec3 pseudo(const quat&) const;
		quat& randu(double = -1, double = +1);

		quat conjugate(void) const;

		vec3 rotate(const vec3&) const;
		vec3 conjugate(const vec3&) const;
		quat conjugate(const quat&) const;

		quat slerp(const quat&, double) const;

		mat3 rotation(void) const;
		double* rotation(double*) const;
	};
}
#pragma once

//std
#include <cstdint>

//Math
#include "Math/inc/Linear/Vector.hpp"

namespace math
{
	class Vec3;
	class Mat3;
}

namespace math
{
	class Quat : public Vector
	{
	public:
		//constructors
		Quat(void);
		Quat(double*);
		Quat(const Quat&);
		Quat(const double*);
		Quat(double, uint32_t);
		Quat(double, const Vec3&);
		Quat(double, double, double, double);
		Quat(const Vec3&, const Vec3&, const Vec3&);
		Quat(const double*, const double*, const double*);

		//destructor
		~Quat(void);

		//operators
		Quat& operator=(const Quat&);
		Quat& operator*=(const Quat&);

		Quat operator*(const Quat&) const;

		double& operator[](uint32_t);
		double& operator()(uint32_t);

		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t) const;

		//views
		Quat& reset(void);
		Quat& view_x(void);
		Quat& view_y(void);
		Quat& view_z(void);
		Quat& isometric(uint32_t);

		//linear
		Quat& normalize(void);
		Vec3 axial(void) const;
		Vec3 pseudo(void) const;
		double angle(void) const;
		Vec3 pseudo(const Vec3&) const;
		Vec3 pseudo(const Quat&) const;
		Quat& randu(double = -1, double = +1);

		Quat conjugate(void) const;

		Vec3 rotate(const Vec3&) const;
		Vec3 conjugate(const Vec3&) const;
		Quat conjugate(const Quat&) const;

		Quat slerp(const Quat&, double) const;

		Mat3 rotation(void) const;
		double* rotation(double*) const;
	};
}
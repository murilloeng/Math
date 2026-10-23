#pragma once

//std
#include <cstdint>

//Math
#include "Math/inc/Linear/Vector.hpp"

namespace math
{
	template <class Group, class Algebra> class LieAlgebra
	{
	public:
		//constructor
		LieAlgebra(void) = default;

		//destructor
		virtual ~LieAlgebra(void) = default;

		//group
		virtual Group exp(void) const = 0;
		virtual Algebra bracket(const Algebra&) const = 0;

		//random
		virtual Algebra& randu(double, double) = 0;

		//linear
		virtual Vector vector(void) = 0;
		virtual Vector vector(void) const = 0;
		virtual Matrix matrix(void) const = 0;
		virtual Matrix adjoint(void) const = 0;
		virtual Matrix jacobian(void) const = 0;
		virtual Matrix jacobian_inverse(void) const = 0;
		virtual Matrix jacobian_gradient(const Vector&) const = 0;

		//operators
		virtual Algebra& operator*=(double) = 0;
		virtual Algebra& operator/=(double) = 0;
		virtual Algebra& operator+=(const Algebra&) = 0;
		virtual Algebra& operator-=(const Algebra&) = 0;

		virtual Algebra operator+(void) const = 0;
		virtual Algebra operator-(void) const = 0;
		virtual Algebra operator*(double) const = 0;
		virtual Algebra operator/(double) const = 0;
		virtual Algebra operator+(const Algebra&) const = 0;
		virtual Algebra operator-(const Algebra&) const = 0;

		virtual double& operator[](uint32_t) = 0;
		virtual const double& operator[](uint32_t) const = 0;
	};
}
#pragma once

//std
#include <cstdint>

//Math
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	template <class Group, class Algebra> class LieAlgebra
	{
	public:
		//constructor
		LieAlgebra(void) = default;

		//destructor
		virtual ~LieAlgebra(void) = default;

		//data
		virtual double* data(void) = 0;
		virtual const double* data(void) const = 0;
		virtual uint32_t dimension(void) const = 0;

		//group
		virtual Group exp(void) const = 0;
		virtual Matrix matrix(void) const = 0;
		virtual Matrix adjoint(void) const = 0;
		virtual Matrix gradient(void) const = 0;
		virtual Matrix gradient_inverse(void) const = 0;
		virtual Algebra bracket(const Algebra&) const = 0;

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
	};
}
#pragma once

//std
#include <cstdint>

//Math
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	template <class Group, class Algebra, uint32_t n, uint32_t d> class LieAlgebra
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
		virtual Vector<d> vector(void) const = 0;
		virtual Matrix<n, n> matrix(void) const = 0;
		virtual Matrix<d, d> adjoint(void) const = 0;
		virtual Matrix<d, d> jacobian(void) const = 0;
		virtual Matrix<d, d> jacobian_inverse(void) const = 0;
		virtual Matrix<d, d> jacobian_gradient(const Vector<d>&) const = 0;

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
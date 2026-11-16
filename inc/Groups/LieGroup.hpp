#pragma once

//std
#include <cstdint>

//Math
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	template <class Group, class Algebra, uint32_t n, uint32_t k> class LieGroup
	{
	public:
		//constructor
		LieGroup(void) = default;

		//destructor
		virtual ~LieGroup(void) = default;

		//random
		virtual Group& randu(double, double) = 0;

		//actions
		virtual Matrix<n, n> matrix(void) const = 0;
		virtual Vector<k> action(const Vector<k>&) const = 0;

		//operations
		virtual Algebra log(void) const = 0;
		virtual Group inverse(void) const = 0;
		virtual Group operator*(const Group&) const = 0;
		virtual Group interpolate(const Group&, double) const = 0;

		virtual double& operator[](uint32_t) = 0;
		virtual const double& operator[](uint32_t) const = 0;
	};
}
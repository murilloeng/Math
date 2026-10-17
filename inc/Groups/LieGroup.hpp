#pragma once

//Math
#include "Math/inc/Linear/Vector.hpp"

namespace math
{
	template <class Group, class Algebra> class LieGroup
	{
	public:
		//constructor
		LieGroup(void) = default;

		//destructor
		virtual ~LieGroup(void) = default;

		//actions
		virtual Matrix matrix(void) const = 0;
		virtual Vector action(const Vector&) const = 0;

		//operations
		virtual Algebra log(void) const = 0;
		virtual Group inverse(void) const = 0;
		virtual Group identity(void) const = 0;
		virtual Group operator*(const Group&) const = 0;
		virtual Group interpolate(const Group&, double) const = 0;
	};
}
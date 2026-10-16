#pragma once

//Math
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	template <class Group, class Algebra> class LieGroup
	{
	public:
		//constructor
		LieGroup(void) = default;

		//destructor
		virtual ~LieGroup(void) = default;

		//operations
		virtual Algebra log(void) const = 0;
		virtual Matrix matrix(void) const = 0;
		virtual Group inverse(void) const = 0;
		virtual Group identity(void) const = 0;
		virtual Group operator*(const Group&) const = 0;
		virtual void action(double*, const double*) const = 0;
		virtual Group interpolate(const Group&, double) const = 0;
	};
}
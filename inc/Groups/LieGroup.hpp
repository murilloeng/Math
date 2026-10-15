#pragma once

//std
#include <memory>

namespace math
{
	class LieAlgebra;
}

namespace math
{
	class LieGroup
	{
	public:
		//constructor
		LieGroup(void);

		//destructor
		~LieGroup(void);

		//operations
		virtual std::unique_ptr<LieAlgebra> log(void) const = 0;
		virtual std::unique_ptr<LieGroup> inverse(void) const = 0;
		virtual std::unique_ptr<LieGroup> identity(void) const = 0;
		virtual std::unique_ptr<LieGroup> operator*(const LieGroup&) const = 0;

		//interpolation
		void interpolate(LieGroup&, const LieGroup&, double t);
	};
}
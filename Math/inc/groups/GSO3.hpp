#pragma once

//math
#include "Math/Math/inc/linear/quat.hpp"

namespace math
{
	namespace groups
	{
		class ASO3;
	}
}

namespace math
{
	namespace groups
	{
		class GSO3
		{
		public:
			//constructor
			GSO3(void);

			//destructor
			~GSO3(void);

			//logarithm
			ASO3 logarithm(void) const;

		private:
			//data
			quat m_quaternion;
		};
	}
}
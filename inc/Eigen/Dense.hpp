#pragma once

//Math
#include "Math/inc/Eigen/Base.hpp"

namespace math
{
	namespace eigen
	{
		class Dense : public Base
		{
		public:
			//constructor
			Dense(uint32_t);

			//destructor
			virtual ~Dense(void);
		};
	}
}
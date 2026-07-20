#pragma once

//Math
#include "Math/inc/Eigen/Eigen.hpp"

namespace math
{
	namespace eigen
	{
		class Dense : public Eigen
		{
		public:
			//constructor
			Dense(uint32_t);

			//destructor
			virtual ~Dense(void);
		};
	}
}
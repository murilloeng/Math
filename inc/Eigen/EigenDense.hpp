#pragma once

//Math
#include "Math/inc/Eigen/Eigen.hpp"

namespace math
{
	namespace eigen
	{
		class EigenDense : public Eigen
		{
		public:
			//constructor
			EigenDense(uint32_t);

			//destructor
			virtual ~EigenDense(void);
		};
	}
}
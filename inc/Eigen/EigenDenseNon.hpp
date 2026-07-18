#pragma once

//Math
#include "Math/inc/Eigen/EigenDense.hpp"

namespace math
{
	namespace eigen
	{
		class EigenDenseNon : public EigenDense
		{
		public:
			//constructor
			EigenDenseNon(uint32_t);

			//destructor
			virtual ~EigenDenseNon(void);
		};
	}
}
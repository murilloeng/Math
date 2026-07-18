#pragma once

//Math
#include "Math/inc/Eigen/Dense.hpp"

namespace math
{
	namespace eigen
	{
		class DenseNon : public Dense
		{
		public:
			//constructor
			DenseNon(uint32_t);

			//destructor
			virtual ~DenseNon(void);
		};
	}
}
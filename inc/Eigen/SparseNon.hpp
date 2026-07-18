#pragma once

//Math
#include "Math/inc/Eigen/Sparse.hpp"

namespace math
{
	namespace eigen
	{
		class SparseNon : public Sparse
		{
		public:
			//constructor
			SparseNon(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*);

			//destructor
			virtual ~SparseNon(void);
		};
	}
}
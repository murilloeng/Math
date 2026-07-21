#pragma once

//Math
#include "Math/inc/Eigen/Sparse.hpp"

namespace math
{
	namespace eigen
	{
		class SparseSym : public Sparse
		{
		public:
			//constructor
			SparseSym(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*);

			//destructor
			virtual ~SparseSym(void);

			//data
			double shift(double);
			double shift(void) const;

		protected:
			//data
			double m_shift;
		};
	}
}
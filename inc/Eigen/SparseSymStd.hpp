#pragma once

//Math
#include "Math/inc/Eigen/SparseSym.hpp"

namespace math
{
	namespace eigen
	{
		class SparseSymStd : public SparseSym
		{
		public:
			//constructor
			SparseSymStd(uint32_t, uint32_t, uint32_t, double*, double*, double*);

			//destructor
			virtual ~SparseSymStd(void);

		protected:
			//data
			double* m_A;
			double* m_s;
			double* m_U;
		};
	}
}
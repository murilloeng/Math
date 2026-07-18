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
			SparseSymStd(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*, const double*, double*, double*);

			//destructor
			virtual ~SparseSymStd(void);

			//compute
			bool compute(void) override;

		protected:
			//operation
			void sparse_product(double*, const double*) const;

			//data
			double* m_s;
			double* m_U;
			const double* m_A;
		};
	}
}
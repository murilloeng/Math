#pragma once

//Math
#include "Math/inc/Eigen/DenseSym.hpp"

namespace math
{
	namespace eigen
	{
		class DenseSymStd : public DenseSym
		{
		public:
			//constructor
			DenseSymStd(uint32_t, double*, double*, double*);
			DenseSymStd(uint32_t, double*, double*, double*, double, double);
			DenseSymStd(uint32_t, double*, double*, double*, uint32_t, uint32_t);

			//destructor
			virtual ~DenseSymStd(void);

			//compute
			bool compute(void) override;

		protected:
			//compute
			bool compute_full(void);
			bool compute_partial(void);

			//data
			double* m_A;
			double* m_s;
			double* m_U;
		};
	}
}
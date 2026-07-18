#pragma once

//Math
#include "Math/inc/Eigen/EigenDenseSym.hpp"

namespace math
{
	namespace eigen
	{
		class EigenDenseSymStd : public EigenDenseSym
		{
		public:
			//constructor
			EigenDenseSymStd(uint32_t, double*, double*, double*);
			EigenDenseSymStd(uint32_t, double*, double*, double*, double, double);
			EigenDenseSymStd(uint32_t, double*, double*, double*, uint32_t, uint32_t);

			//destructor
			virtual ~EigenDenseSymStd(void);

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
#pragma once

//Math
#include "Math/inc/Eigen/EigenDenseSym.hpp"

namespace math
{
	namespace eigen
	{
		class EigenDenseSymGen : public EigenDenseSym
		{
		public:
			//constructor
			EigenDenseSymGen(uint32_t, double*, double*, double*, double*);

			//destructor
			virtual ~EigenDenseSymGen(void);

			//compute
			bool compute(void) override;

		protected:
			//compute
			bool compute_full(void);
			bool compute_partial(void);

			//data
			double* m_A;
			double* m_B;
			double* m_s;
			double* m_U;
		};
	}
}
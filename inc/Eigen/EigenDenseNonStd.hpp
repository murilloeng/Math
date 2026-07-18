#pragma once

//Math
#include "Math/inc/Eigen/EigenDenseNon.hpp"

namespace math
{
	namespace eigen
	{
		class EigenDenseNonStd : public EigenDenseNon
		{
		public:
			//constructor
			EigenDenseNonStd(uint32_t, double*, double*, double*, double*, double*);

			//destructor
			virtual ~EigenDenseNonStd(void);

			//compute
			bool compute(void) override;

		protected:
			//data
			double* m_A;
			double* m_U;
			double* m_V;
			double* m_sr;
			double* m_si;
		};
	}
}
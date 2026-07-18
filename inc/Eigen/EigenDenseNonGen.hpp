#pragma once

//Math
#include "Math/inc/Eigen/EigenDenseNon.hpp"

namespace math
{
	namespace eigen
	{
		class EigenDenseNonGen : public EigenDenseNon
		{
		public:
			//constructor
			EigenDenseNonGen(uint32_t, double*, double*, double*, double*, double*, double*);

			//destructor
			virtual ~EigenDenseNonGen(void);

			//compute
			bool compute(void) override;

		protected:
			//data
			double* m_A;
			double* m_B;
			double* m_U;
			double* m_V;
			double* m_sr;
			double* m_si;
		};
	}
}
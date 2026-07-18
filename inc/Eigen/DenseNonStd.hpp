#pragma once

//Math
#include "Math/inc/Eigen/DenseNon.hpp"

namespace math
{
	namespace eigen
	{
		class DenseNonStd : public DenseNon
		{
		public:
			//constructor
			DenseNonStd(uint32_t, double*, double*, double*, double*, double*);

			//destructor
			virtual ~DenseNonStd(void);

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
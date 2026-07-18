#pragma once

//Math
#include "Math/inc/Eigen/DenseNon.hpp"

namespace math
{
	namespace eigen
	{
		class DenseNonGen : public DenseNon
		{
		public:
			//constructor
			DenseNonGen(uint32_t, double*, double*, double*, double*, double*, double*);

			//destructor
			virtual ~DenseNonGen(void);

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
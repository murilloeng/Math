#pragma once

//Math
#include "Math/inc/Eigen/EigenDense.hpp"

namespace math
{
	namespace eigen
	{
		class EigenDenseStd : public EigenDense
		{
		public:
			//constructor
			EigenDenseStd(uint32_t, double*);

			//destructor
			virtual ~EigenDenseStd(void);

		protected:
			//data
			double* m_A;
		};
	}
}